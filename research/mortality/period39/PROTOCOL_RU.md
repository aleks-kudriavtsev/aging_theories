#39: российский возрастно-половой периодный ориентир

Протокол фиксируется до новых5/10-летних численных результатов. Источник — уже сохранённый и hash-проверенный derivedCSV UNWPP2024, Россия2023, годовые возраста40–89, полF/M. Это не новые individual outcomes и не зарегистрированные смерти2023. OriginalRDA checksum и числовое scaling заново не аутентифицированы; используем явно названный local-derived reference по документированным единицам.

Для исходного точного возрастаa на день рождения и горизонтаh берём возрастаa…a+h−1. Поэтому старт80/10лет требует80…89, endpoint90. Фиксируем schedule2023 на все будущие возрастные годы: это условный period reference, не прогноз реальных календарных2026–2036 и не персональный клинический риск.

Primaryconstantforce:Q=1−exp(−Σmx). SecondaryUDD a_x=.5:qx=mx/(1+.5mx),Q=1−Π(1−qx). Это разные допущения, не статистически отобранные модели и не взаимозаменяемые формулы. Сохранённый авторскийqx — отдельный diagnostic, не официальный WPPqx.

Заранее задаём164age/sex/horizon cells,20selectedreportcells,72chaining checks. Два расчёта чувствительности к округлениюmx ±.5e−6 — арифметические bounds при этом допущении, не sampling/forecastCI; отдельная model-assumption difference CFvsUDD не confidence interval. Не вводим smoothage fits и не выбираем вариант по несуществующим validation outcomes.

Строгий adapter требует exactbirthday, ageinteger40–80, sexF/M, RUS,source2023,horizon5/10,явное period assumption и purposepopulation_reference_only. Он не принимает biomarker/disease-free/stagebaseline запрос и не игнорирует неизвестные поля. Общая популяционная смертность уже включает больных; её нельзя прибавлять к postdiagnosis allcause curves или выдавать за diseasefree intensity. Для совместной calibration нужны representative riskfactor/state distribution и независимые outcomes.

Параллельная SHFM qualification — только primary source/rights/version review. Численная5-летняя SHFM reproduction не выполняется без отдельного completecontract и frozenprotocol. Все прежние модели/результаты38ираньше сохраняются; worklocalonly,no push/PR/merge/deploy.
