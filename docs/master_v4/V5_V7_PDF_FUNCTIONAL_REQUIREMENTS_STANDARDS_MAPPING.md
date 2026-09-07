# تحليل `نماذج V5-V7.pdf` — المتطلبات الوظيفية ومطابقتها بالمعايير الحالية

## 1. النتيجة العامة

ملف `نماذجV5-V7.pdf` هو **ملف لوحات هندسية مرئية image-only** مكوّن من 15 صفحة. لم ينتج `pdftotext` أي نص قابل للاستخراج، وتبيّن أن كل صفحة تقريباً صورة JPEG مضمّنة. لذلك تم تحليل اللوحات بصرياً، مع الحفاظ على تصنيفها كـ**REFERENCE / DESIGN INTENT**.

تتضمن اللوحات رسومات exploded views، مخططات مسارات حمل، جداول مواد وأجزاء، أبعاداً مرئية، force–stroke curves، حالات lock/unlock، مراحل زمنية للحادث، ومخططات تحقق. لكنها لا تقدم وحدها native CAD، أو PMI/GD&T قابلاً للقراءة الآلية، أو شهادات مواد، أو بيانات معايرة، أو raw test data، أو تقرير مطابقة تنظيمي.

**القرار:** توفر اللوحات أساساً جيداً لاستخراج متطلبات الهندسة الوظيفية وبناء خطة تحقق، لكنها لا تغلق أي متطلب تصنيع أو اعتماد. الحالة الحالية هي:

```text
FUNCTIONAL REQUIREMENTS EXTRACTED
STANDARDS ALIGNMENT PARTIAL
COMPLIANCE NOT ESTABLISHED
STOP RELEASE
```

## 2. منهج التحليل

تم فحص metadata وembedded images لكل الصفحات. كما تمت مراجعة اللوحات صفحةً صفحة. تم تسجيل كل متطلب بصيغة:

`الوظيفة → متغير قابل للقياس → دليل مطلوب → معيار محتمل → حالة الدليل الحالية`.

الأرقام الظاهرة على اللوحات، مثل `1120 mm` و`560 mm` و`520 mm` و`460 mm`، وكذلك الشوط أو منحنيات القوة، تعاملت معها كـ**قيم مرجعية مرئية** فقط. لم يتم تحويلها إلى nominal release أو acceptance criterion.

## 3. المتطلبات المستخرجة حسب الصفحات

| المعرّف | الصفحات | ما يظهر في اللوحة | المتطلب الوظيفي المستخلص | الدليل اللازم |
|---|---:|---|---|---|
| V5F-01 | 1–5 | هيكل مقعد، rail pair، seatback، absorber، pan، restraint | تكامل مقعد مع مساري حمل طوليين ومسار امتصاص وطريق دوران منفصل | native assembly، joints، hardpoints، interference، seat/rail load test |
| V5F-02 | 6 | exploded lock/unlock mechanism، locked/unlocked sections، trigger | القفل يجب أن ينتقل بين حالات محددة، ويظل مقفلاً تحت الحمل، ولا يحرر دون trigger مقصود | state-cycle، load-hold، unintended-release، reset، failure-mode test |
| V5F-03 | 7 | absorber cartridge، force–stroke curve، stroke cushion | وحدة امتصاص قابلة للاستبدال تنتج قانون قوة–شوط قابلاً للقياس ومتكرراً | calibrated force/displacement، hysteresis، temperature، cycling، serial traceability |
| V5F-04 | 7 | rail system وconnection details | توجيه المقعد على rail pair مع منع الانحراف أو الانفصال وتوفير end stops | parallelism، straightness، retention، stop-load، left/right reaction measurement |
| V5F-05 | 7 | welding specification وmaterial/performance tables | الوصلات واللحامات يجب أن تنقل الحمل المطلوب دون فشل غير مقبول | released weld specification، WPS/PQR عند الحاجة، material certificates، inspection criteria |
| V5F-06 | 8 | normal → detection → pelvis control → energy management → rebound → post-event | نظام زمني متعدد المراحل يغيّر القيود الميكانيكية وفق حالة الحدث | synchronized acceleration/travel/angle/state channels، trigger latency، phase acceptance limits |
| V5F-07 | 8 | force/time curves وload-path overview | توزيع الحمل بين المسارات يجب أن يبقى قابلاً للقياس، مع تحديد contributions لكل مسار | channel-level calibration، CFC/filter، timebase، uncertainty، raw-data hash |
| V5F-08 | 9 | front/side/rear views، sections، dimensions، general tolerances، belt/seat interfaces | وجود envelope وdatums وinterfaces وتفاوتات تصنيع قابلة للفحص | released drawing، native CAD، PMI/GD&T، tolerance stack، inspection report |
| V5F-09 | 10 | detailed exploded CAD board، joints/details، material specifications | assembly decomposition مع joints قابلة للخدمة ووحدات قابلة للاستبدال | BOM، mating conditions، fastener identity/preload، service access، revision control |
| V5F-10 | 11–12 | rear-seat configurations، industrial design، upholstery وfinish | تكامل المقعد داخل rear-seat packaging مع الحفاظ على الوظائف الهيكلية والتقييدية | vehicle hardpoint definition، occupant clearance، belt routing، trim interference، regulatory test configuration |
| V5F-11 | 13 | أبعاد مرئية تقريباً 1120/560/520/460 mm، phase timing bands، state table، force–stroke | تعريف package envelope وشروط الحالة والتوقيت والخرطوشة كمدخلات تصميمية | dimension source/revision، tolerance، state-machine specification، calibrated dynamic validation |
| V5F-12 | 14–15 | multi-domain mechanism، sensor/actuator، control zones، occupant models، validation plan | ربط sensors/actuators بالحالات الميكانيكية والتحكم في trajectory دون ادعاء حماية مثبتة | sensor identity/calibration، actuator fail-safe، trigger repeatability، ATD/sled evidence، post-test inspection |

## 4. تحليل خاص بآلية lock/unlock — الصفحة 6

تُظهر الصفحة 6 آلية قفل متعددة الأجزاء، حالات `locked/unlocked`، آلية trigger، وتسلسلاً للحركة. هندسياً، ينتج ذلك المتطلبات التالية:

| المتطلب | معيار القبول الذي يجب تعريفه | العلاقة بالمعايير |
|---|---|---|
| positive engagement | لا يحدث release تحت الحمل المحدد | FMVSS 207 للمقاعد القابلة للطي/مسند الظهر، وUNECE R17 لأنظمة القفل والضبط بحسب التطبيق |
| intentional release | لا يتحقق التحرير إلا من trigger مقصود وتحت شروطه | لا يحدد FMVSS 207 وحده كل منطق trigger الخاص بـV5؛ يحتاج مواصفة داخلية |
| repeatability | عدد دورات محدد دون زيادة play أو failure | متطلب تحقق داخلي، وليس حكماً تنظيمياً عاماً |
| reset/serviceability | إعادة الضبط يمكن تنفيذها دون فقدان حالة القفل أو التقييد | متطلب هندسي داخلي مع inspection procedure |
| fault containment | spring/pin/trigger failure لا يؤدي إلى release غير متحكم به | FMEA واختبار fault insertion مطلوبان |

FMVSS 207 يتناول self-locking restraint للمقاعد أو مساند الظهر القابلة للطي في الحالات المنطبقة، ويتضمن قوى طولية واختبار تسارع 20 g للجزء القابل للطي [3]. أما UNECE R17 فيتناول seat adjustment and locking systems ضمن مجال المقاعد وتثبيتها [8]. هذا لا يثبت أن آلية V5 تحقق تلك المتطلبات؛ يلزم تعريف المركبة، اتجاه المقعد، موضعه، وحالة الضبط ثم اختبارها.

## 5. تحليل absorber وforce–stroke — الصفحتان 7 و13

تُظهر اللوحات cartridge قابلة للاستبدال ومنحنى force–stroke ومراحل امتصاص. هذه الوظيفة تتطلب فصل أربعة أشياء:

1. **قانون التصميم:** ما القوة المستهدفة وما الشوط المقصود.
2. **قانون الجزء الفعلي:** force–stroke وhysteresis وtemperature وrate dependence.
3. **تفاعل النظام:** أثر الخرطوشة على rail loads وseatback rotation وrestraint path.
4. **قبول النظام:** هل يبقى المقعد والتثبيت ضمن متطلبات FMVSS/UNECE والبروتوكول الديناميكي المعين؟

FMVSS 207 وUNECE R17 قد يفرضان متطلبات قوة المقعد وتثبيته، لكنهما لا يصادقان تلقائياً على absorber law أو على الشوط المرئي في اللوحة. لذلك لا يمكن استخدام منحنى اللوحة كبديل عن characterization test. كما أن قيم `18–22 kN` أو أي منحنى مرئي يجب أن تبقى design reference إلى أن ترتبط بحساس معاير وraw data وconfiguration hash.

## 6. تحليل crash load path والتوقيت — الصفحة 8

تقدم الصفحة 8 تسلسلاً زمنياً تقريبياً من الحالة الطبيعية إلى detection ثم pelvis control وenergy management وrebound وpost-event. هذا يخلق متطلبات state-machine قابلة للتحقق:

| حقل مطلوب | سبب الحاجة |
|---|---|
| الحالة الابتدائية | لمنع تفسير أي trigger على أنه حادث |
| signal/threshold | لتحديد سبب الانتقال وليس مجرد ملاحظة لاحقة |
| latency | لقياس ما إذا كان الانتقال يحدث قبل أو بعد الحركة الحرجة |
| mechanical constraint per state | لإثبات أن الحالة تغيّر الحمل فعلياً |
| travel/rotation limit | لمنع تجاوز envelope |
| unintended release criterion | لإثبات fail-safe |
| sensor time base and synchronization | لربط acceleration وtravel وangle وforce |
| post-event state | لتحديد إمكانية البقاء الآمن أو إعادة الضبط |

ISO 6487 وSAE J211/1 مناسبان لإطار القياس، لأنهما يتعاملان مع قناة القياس كاملة، والتوقيت، والاستجابة الترددية، والمعايرة، والترشيح، وعدم اليقين [11] [12]. لكنهما لا يحددان وحدهما نجاح state-machine أو سلامة occupant.

## 7. تحليل لوحات الرسم والتصنيع — الصفحات 9 و10 و12 و13

تظهر اللوحات views وsections وBOM/material tables وgeneral tolerances وweld/join callouts وseat-belt interfaces. هذه العناصر هي متطلبات **تعريف تصنيع** وليست دليلاً على أن التصنيع أصبح معتمداً.

| عنصر اللوحة | ما يجب أن يكون موجوداً للإفراج الهندسي | الحالة الحالية |
|---|---|---|
| geometry | native CAD وrevision وconfiguration ID | مفقود |
| dimensions | dimension source وdatums وunits وrevision | مرئي بصرياً فقط |
| tolerances | GD&T وtolerance stack وinspection method | غير قابل للتحقق من raster board |
| materials | grade، heat treatment، certificate، lot mapping | غير مثبت |
| welds/joints | WPS/PQR أو joining specification وinspection class | غير مثبت |
| fasteners | standard، grade، preload، locking method، traceability | غير مثبت |
| BOM | released part numbers، quantities، alternates، authority | reference board فقط |
| belt/seat interface | anchorage geometry and strength evidence | غير مثبت |
| service access | removal path، replaceable module، safe reset | concept-level |

لذلك فإن كلمة `production ready` الظاهرة في بعض اللوحات يجب قراءتها كعنوان للوحة أو هدف تصميمي، وليس كإثبات لإصدار تصنيع.

## 8. المطابقة مع المعايير الحالية

### FMVSS 207 — المقاعد وأنظمة تثبيتها

يرتبط مباشرة بوظائف المقعد والـrail والـseatback والـlocking. يفرض، ضمن نطاقه واستثناءاته، قوى طولية أمامية وخلفية مقدارها `20 × كتلة المقعد × 9.8 N`، وعزماً قدره `373 N·m` لكل موضع جلوس، ومتطلبات الاحتفاظ بالضبط. كما يعالج self-locking للمقاعد أو مساند الظهر القابلة للطي في الحالات المنطبقة [3].

**المطابقة:** `PARTIAL`. اللوحات تحدد الوظائف التي يجب اختبارها، لكنها لا تحدد فئة المركبة، حالة التصنيع، configuration المعتمدة، أو نتائج الاختبار.

### FMVSS 210 — تثبيتات أحزمة الأمان

يرتبط مباشرة بالـseat-belt interface الظاهر في الصفحات 9–13 وباستمرار مسار التقييد أثناء حركة المقعد. النص الحالي يتضمن 5,000 lb لبعض حالات التثبيت و3,000 lb لقوى lap/shoulder المتزامنة في حالات أخرى، مع متطلبات موضعية وإجراءات body-block/FAD واختبار hold [4].

**المطابقة:** `INTERFACE RELEVANCE ONLY`. وجود belt interface في اللوحة لا يثبت قوة anchorage أو صلاحية مسار الحزام أثناء ride-down.

### FMVSS 202a — مساند الرأس

يرتبط بالـseatback والـhead/torso guidance في الصفحات 1–5 و10–15. يغطي الارتفاع والعرض والـbackset والفجوات وامتصاص الطاقة والاحتفاظ بالضبط والقوة. ويوفر مساراً ديناميكياً بحد دوران خلفي 12 درجة و`HIC15 ≤ 500`، مع شروط الاختبار الأخرى [5].

**المطابقة:** `PARTIAL`. صور المقعد والـhead support لا تكفي لإثبات أن العنصر مسند رأس مطابق أو أن trajectory control يحقق معيار إصابة.

### UNECE R14 وR17 وR25

| المعيار | الصلة باللوحات | النتيجة |
|---|---|---|
| R14 | belt anchorages، مواقعها وقوتها، بما فيها seat-mounted anchorages عند انطباقها | مطلوب released anchorage package واختبار القوة/الموقع |
| R17 | seat strength، seat anchorages، adjustment/locking، وبعض head-restraint provisions | مطلوب تحديد فئة المركبة واتجاه المقعد configuration ثم تنفيذ الاختبار |
| R25 | head-restraint device عندما لا يكون مسار R17 هو المسار الحاكم | يحتاج تحديد قاعدة عدم الازدواج والنسخة التنظيمية الفعلية |

لا يوجد معيار واحد يغطي كل الوظائف الظاهرة في اللوحات. يجب بناء compliance matrix حسب سوق التطبيق وفئة المركبة وسنة التصنيع وسلسلة التعديلات السارية [6] [7] [8].

### ISO 6487 وSAE J211/1

تصلح هذه المعايير لمتطلبات الصفحتين 8 و13 المتعلقة بالقنوات، force/time، acceleration، sensor/actuator timing، وvalidation plan. وهي تغطي خصائص قناة القياس، المعايرة، CFC، phase delay، timebase، uncertainty، والتقرير [9] [10]. لكنها لا تقرر أن المنحنى أو seat load أو injury metric ناجح؛ معيار القبول يأتي من اللائحة أو بروتوكول الاختبار المحدد.

## 9. القيم المرئية التي لا يجوز ترقيتها إلى متطلبات اعتماد

| القيمة أو المحتوى | تصنيف صحيح |
|---|---|
| `1120 mm`, `560 mm`, `520 mm`, `460 mm` الظاهرة على الصفحة 13 | reference dimensions تحتاج مصدر رسم أصلي وتفاوتات |
| force–stroke curves | design/reference curves حتى توفر raw calibrated data |
| phase timing bands | proposed control strategy حتى توفر state-machine acceptance limits |
| material/color/finish tables | design/material intent حتى توفر certificates وlot mapping |
| “production ready” أو “industrial design” | board label، لا release authority |
| CAD format icons مثل STEP/IGES/CATIA/SolidWorks | إشارة إلى صيغ مرغوبة أو متوقعة، لا دليل على توفر الملفات |
| crash sequence 0–400 ms | scenario/illustration حتى توفر measured synchronized channels |
| occupant models وHIC/Nij references | analysis concept حتى يتم تنفيذ test/CAE correlation |

## 10. مصفوفة الأدلة المطلوبة

| بوابة | دليل مطلوب من اللوحات | معيار/إطار مساند | الحالة الحالية |
|---|---|---|---|
| CAD | native parametric assembly، interference، travel، envelope | project gate / FMVSS 207 preparation | BLOCKED |
| Lock | state cycle، no unintended release، load/acceleration hold | FMVSS 207 / UNECE R17 where applicable | BLOCKED |
| Absorber | calibrated force–stroke، hysteresis، temperature، cycling | project-specific; related seat-strength tests | BLOCKED |
| Anchoring | location، hardware، simultaneous lap/shoulder، strength | FMVSS 210 / UNECE R14 | BLOCKED |
| Head restraint | height، width، backset، gaps، strength، energy، dynamic | FMVSS 202a / UNECE R17/R25 | BLOCKED |
| Instrumentation | calibration، CFC، timebase، sync، uncertainty، raw hash | ISO 6487 / SAE J211/1 | CONTRACT ONLY |
| Dynamic validation | ATD/sled kinematics، loads، post-test inspection | applicable regulation/protocol | NOT AVAILABLE |
| Release | revision، BOM، materials، PMI/GD&T، inspection | manufacturer quality system | NOT AUTHORIZED |

## 11. الاستنتاج

يستخرج `نماذجV5-V7.pdf` متطلبات وظيفية مفيدة، خصوصاً في مجالات lock/unlock، absorber، rail guidance، time-programmed state transitions، seat/vehicle integration، sensor/actuator control، والـvalidation planning. لكنه يبقى **مرجع تصميم مرئي**.

المطابقة الحالية مع FMVSS وUNECE وISO/SAE هي مطابقة موضوعية ومتطلباتية، وليست إعلان compliance. لا يمكن إغلاق المتطلبات دون native CAD، BOM مفرج عنه، مواد وشهادات، PMI/GD&T، معايرة، raw data، اختبارات قوة وتثبيت ومساند رأس، ونتائج sled/CAE قابلة لإعادة الإنتاج.

يبقى القرار: **STOP RELEASE**. لا يجوز إنشاء manufacturer-release package أو تحويل اللوحات إلى تصنيع أو اعتماد قبل استرداد المصدر الأصلي المفرج عنه وإعادة بناء evidence graph.

## References

[1]: file:///home/ubuntu/upload/نماذجV5-V7.pdf "T-OCS V5–V7 Engineering Forms and Visual Boards"

[2]: file:///home/ubuntu/upload/T_OCS_V7_Patent_Draft_Invention_Disclosure_.pdf "T-OCS V7.0 Patent Preparation — V7 Enhancement Addendum"

[3]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.207 "49 CFR 571.207 — Standard No. 207; Seating systems"

[4]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.210 "49 CFR 571.210 — Standard No. 210; Seat belt assembly anchorages"

[5]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.202a "49 CFR 571.202a — Standard No. 202a; Head restraints"

[6]: https://unece.org/sites/default/files/2025-05/R014r6e%20%282%29.pdf "UNECE UN Regulation No. 14 — Safety-belt anchorages"

[7]: https://unece.org/sites/default/files/2022-01/R017r6e.pdf "UNECE UN Regulation No. 17 — Seats, their anchorages and head restraints"

[8]: https://unece.org/sites/default/files/2025-08/R025r1a4e.pdf "UNECE UN Regulation No. 25 — Head restraints"

[9]: https://www.iso.org/standard/64041.html "ISO 6487:2015 — Road vehicles, measurement techniques in impact tests, instrumentation"

[10]: https://saemobilus.sae.org/standards/j2111_202208-instrumentation-impact-test-part-1-electronic-instrumentation "SAE J211/1_202208 — Instrumentation for Impact Test, Part 1: Electronic Instrumentation"
