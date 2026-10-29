"""Synthetic FWF contract tests, without source downloads or participant data."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from research.mortality.routine29.reader import read_checked

FIELDS = {'SEQN': [1,3], 'CREAT': [5,8], 'BMI': [10,13], 'FLAG': [15,15]}
NOMINAL = 20


def dictionary(missover=True, **spans):
    result = "/* metadata only */\nLENGTH\n SEQN 7\n CREAT 8\n BMI 8\n;\n"
    result += 'FORMAT\n CREAT Z6.1\n BMI Z6.1\n;\n'
    result += 'INFILE DATA' + (' MISSOVER' if missover else '') + ' LRECL=22;\nINPUT\n'
    for key,(start,end) in FIELDS.items():
        span=str(start) if start==end else f'{start}-{end}'
        result += f' {key} {spans.get(key, span)}\n'
    return result+';\nLABEL SEQN="Identifier";\n'


def row(**updates):
    values = dict(SEQN='1',CREAT='.4',BMI='27.3',FLAG='0')
    values.update(updates)
    line=[' ']*NOMINAL
    for key,(start,end) in FIELDS.items():
        token=str(values[key]).rjust(end-start+1)
        if len(token)!=end-start+1:
            raise AssertionError('Synthetic test literal is too wide')
        line[start-1:end]=token
    return ''.join(line)


class VerifiedFWFTests(unittest.TestCase):
    def read(self, data, sas=None, expected_rows=None, fields=None, nominal_width=NOMINAL,
             source_newline='\n',sas_newline='\n'):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            content=data.replace('\n',source_newline).encode('ascii')
            (root/'source.dat').write_bytes(content)
            (root/'source.sas').write_bytes((dictionary() if sas is None else sas).replace('\n',sas_newline).encode('ascii'))
            frame,audit=read_checked(root/'source.dat',root/'source.sas',FIELDS if fields is None else fields,
                                     nominal_width,expected_rows)
            return frame,audit,content

    def test_numeric_literals_preserved_without_format_scaling(self):
        d,audit,_=self.read(row()+'\n'+row(SEQN=2,CREAT='-0.1',BMI=8888)+'\n',expected_rows=2)
        np.testing.assert_array_equal(d.CREAT.to_numpy(),[.4,-.1])
        np.testing.assert_array_equal(d.BMI.to_numpy(),[27.3,8888.])
        self.assertFalse(audit['FORMAT_scaling_applied'])
        self.assertEqual(audit['cell_comparisons'],8)
        self.assertTrue(audit['exact_match'])

    def test_only_explicit_dot_and_blank_are_missing(self):
        d,_,_=self.read(row(CREAT='.',BMI='')+'\n'+row(SEQN=2,CREAT='0',BMI='1.')+'\n')
        self.assertTrue(pd.isna(d.CREAT.iloc[0]));self.assertTrue(pd.isna(d.BMI.iloc[0]))
        self.assertEqual(d.CREAT.iloc[1],0.)
        self.assertEqual(d.BMI.iloc[1],1.)

    def test_regular_signed_and_exponent_tokens_are_numbers(self):
        d,_,_=self.read(row(CREAT='+.4',BMI='1e2'))
        self.assertEqual(d.CREAT.iloc[0],.4);self.assertEqual(d.BMI.iloc[0],100.)

    def test_nan_inf_control_characters_and_junk_rejected(self):
        for value in ('NaN','inf','-inf','1e99','1 2','x','.A','\t.4'):
            # 1e99 is finite and legitimate, unlike overflowing exponent.
            if value=='1e99':
                d,_,_=self.read(row(CREAT=value));self.assertEqual(d.CREAT.iloc[0],1e99)
                continue
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError,'source field'):
                    self.read(row(CREAT=value))

    def test_overflow_token_is_rejected(self):
        fields={'SEQN':[1,3],'VALUE':[5,10]}
        sas='INFILE DATA MISSOVER;\nINPUT\n SEQN 1-3\n VALUE 5-10\n;'
        content='  1 1e9999'+' '*10
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            self.read(content,sas=sas,fields=fields)

    def test_input_isolated_from_length_and_format_blocks(self):
        d,_,_=self.read(row(),sas=dictionary(),sas_newline='\r\n')
        self.assertEqual(d.SEQN.iloc[0],1.)
        with self.assertRaisesRegex(ValueError,'Exact INPUT positions'):
            self.read(row(),sas=dictionary(CREAT='6-9'))

    def test_spaced_wrong_range_not_accepted_as_single_column(self):
        with self.assertRaisesRegex(ValueError,'Exact INPUT positions'):
            self.read(row(),sas=dictionary(FLAG='15 - 16'))

    def test_extra_or_duplicate_input_declaration_rejected(self):
        for sas in (dictionary()+'\nINPUT\n SEQN 1-3\n;',
                    dictionary().replace(' CREAT 5-8',' CREAT 5-8\n CREAT 5-8')):
            with self.assertRaises(ValueError):self.read(row(),sas=sas)

    def test_trailing_absence_only_with_missover(self):
        short=row()[:15]
        d,audit,_=self.read(row()+'\n'+row(SEQN=2)[:19]+'\n'+row(SEQN=3)[:15])
        self.assertEqual(len(d),3)
        self.assertEqual(audit['source_record_width_counts'],{15:1,19:1,20:1})
        self.assertEqual(audit['rows_with_absent_trailing_unselected_characters'],2)
        self.assertEqual(audit['selected_fields_end'],15)
        with self.assertRaisesRegex(ValueError,'record width'):
            self.read(short,sas=dictionary(missover=False))

    def test_no_arbitrary_short_overlong_blank_or_eof_record_skip(self):
        for content in (row()[:14],row()+' ',row()+'\n\n',row()+'\n\x1a'):
            with self.subTest(length=len(content)):
                with self.assertRaisesRegex(ValueError,'record width'):self.read(content)

    def test_missover_in_label_or_quoted_filename_does_not_enable_it(self):
        for sas in (dictionary(missover=False)+'\nLABEL BMI="MISSOVER";',
                    dictionary(missover=False).replace('INFILE DATA','INFILE "MISSOVER"')):
            with self.assertRaisesRegex(ValueError,'record width'):
                self.read(row()[:15],sas=sas)

    def test_comment_input_examples_do_not_override_actual_input(self):
        sas='/*\nINPUT\n SEQN 7\n;\n*/\n'+dictionary()
        d,_,_=self.read(row(),sas=sas);self.assertEqual(d.SEQN.iloc[0],1.)

    def test_exact_source_count_is_required_when_supplied(self):
        with self.assertRaisesRegex(ValueError,'locked expected count'):
            self.read(row()+'\n'+row(SEQN=2),expected_rows=3)
        d,audit,_=self.read(row(),expected_rows=1)
        self.assertEqual(len(d),1);self.assertEqual(audit['expected_rows'],1)

    def test_identifiers_missing_duplicate_fractional_or_nonpositive_fail(self):
        for data in (row(SEQN='.'),row()+'\n'+row(),row(SEQN='0'),row(SEQN='1.5')):
            with self.assertRaisesRegex(ValueError,'SEQN'):self.read(data)

    def test_source_hash_counts_original_lf_crlf_and_no_terminal_newline(self):
        for newline in ('\n','\r\n'):
            for final in ('','\n'):
                d,audit,raw=self.read(row()+'\n'+row(SEQN=2)+final,source_newline=newline)
                self.assertEqual(len(d),2)
                self.assertEqual(audit['source_sha256'],hashlib.sha256(raw).hexdigest())
                self.assertEqual(audit['source_bytes'],len(raw))

    def test_independent_reader_mismatch_is_not_silenced(self):
        with patch('research.mortality.routine29.reader.pd.read_fwf',
                   return_value=pd.DataFrame({'SEQN':[1.],'CREAT':[9.],'BMI':[27.3],'FLAG':[0.]})):
            with self.assertRaisesRegex(ValueError,'Independent source readers disagree'):self.read(row())

    def test_bad_field_contracts_and_counts_fail(self):
        for kwargs in ({'fields':{}},{'fields':{'SEQN':[1,3],'BAD':[3,4]}},
                       {'fields':{'SEQN':[False,3]}},{'nominal_width':True},
                       {'fields':{'SEQN':[1,21]}},{'expected_rows':True},
                       {'expected_rows':1.0},{'expected_rows':0}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):self.read(row(),**kwargs)

    def test_empty_file_is_not_valid_empty_release(self):
        with self.assertRaisesRegex(ValueError,'Empty source'):self.read('')

if __name__=='__main__':unittest.main()
