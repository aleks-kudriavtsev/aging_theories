# Source identity for final training (before first final fit)

The CDC index.htm is a live release catalogue, not a participant measurement
file. Earlier cohort27/binary28 gates already excluded this live HTML from
exact historical content equality while checking its official host and resolved
release links. Routine29 validates actual source hashes and declared URLs.

Iteration30 requires exact agreement with iteration29 for **all 13 non-index
files**, including participant data, SAS definitions and source PDF/HTML
documentation. No measurement or codebook change is silently accepted. The live
index is retrieved through the existing official-only route and must resolve
the same release URLs as the saved routine29 SOURCE_MANIFEST. Record both old
and actual index hashes, and assert that it is the only possible changed file.
If any non-index byte hash or release URL changes, stop before final training.

This policy allows catalogue rendering changes, not a different study release.
The actual 14-file retrieval manifest is bound into the new deterministic freeze
manifest. It is fixed before the first final fit and does not depend on outcomes.
