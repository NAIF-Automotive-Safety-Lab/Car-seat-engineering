# تحليل المتطلبات الوظيفية في ملفات PDF ومطابقتها بالمعايير الحالية

## 1. الخلاصة التنفيذية

تقدم ملفات PDF المرفقة **معمارية وظيفية ومركزاً ابتكارياً مقترحاً** لنظام مقعد يتحكم في مسار حركة المقعد، ودوران مسند الظهر، والارتداد، وحركة الحوض، مع الحفاظ على مسار التقييد. لكنها لا تمثل حزمة CAD أصلية مفرجاً عنها، ولا رسومات تصنيع معتمدة، ولا شهادة مطابقة تنظيمية، ولا نتائج اختبار مقاسة.

المعيار الأهم في القراءة الصحيحة للملفات هو أن قيم `25.6 g` و`180 mm` و`18–22 kN` ونتائج `HIC/Nij` الواردة في المواد المصدرية **ليست نتائج T-OCS مقاسة**. ينص ملف V7 نفسه على أنها مدخلات دراسة أو أهداف تصميمية إلى أن يتم التحقق منها باختبار أو تحليل قابل لإعادة الإنتاج.

تتوافق المعمارية الوظيفية جزئياً مع موضوعات المعايير الحالية، لكنها لا تثبت المطابقة. توجد مواءمة قوية مع متطلبات قوة المقاعد وتثبيتاتها في **FMVSS 207**، ومتطلبات تثبيت أحزمة الأمان في **FMVSS 210**، ومتطلبات مساند الرأس في **FMVSS 202a**، ومع تقسيم الوظائف بين **UNECE R14 وR17 وR25**. أما **ISO 6487 وSAE J211/1** فتدعمان متطلبات القياس والمعايرة والتوقيت والترشيح، ولا تحددان وحدهما نجاح أو فشل أداء المقعد.

**الحكم الحالي:** `PARTIAL FUNCTIONAL ALIGNMENT — NOT A COMPLIANCE DETERMINATION`. وتبقى بوابة المشروع الحالية **STOP RELEASE** بسبب غياب حزمة V7 الأصلية وبيانات الاختبار والمواد والتفاوتات والـCAE الجاهز.

## 2. أساس التحليل

تمت قراءة ملف `T_OCS_V7_Patent_Draft_Invention_Disclosure_.pdf`، وبالأخص صفحات V7-1 إلى V7-8، مع استخدام ملف `نماذجV5-V7.pdf` كمرجع بصري للمخططات والأهداف التصميمية. كما تمت مراجعة حالة المستودع وتقارير الاسترداد السابقة، التي تثبت عدم توفر V7 native CAD authoritative package.

> «The strongest V7 independent claim is intentionally narrower than a generic energy attenuating seat claim.» — T-OCS V7، صفحة 3.

هذه العبارة مهمة لأنها تحدد أن التحليل يجب أن يختبر **العلاقة الوظيفية بين المسارات والحالات**، لا أن يعامل كل rail أو absorber أو lock كاختراع مستقل أو كمنتج مطابق.

## 3. المتطلبات الوظيفية المستخرجة من PDF

| المعرّف | المتطلب الوظيفي المستخرج | موضع المصدر | نوع المتطلب | معيار التحقق المقترح |
|---|---|---|---|---|
| FR-01 | وجود مساري حمل طوليين متباعدين جانبياً `110L/110R` لتوجيه الحمل وتقاسمه والتحكم الالتوائي | V7 صفحة 2، V7 صفحة 3 | معماري/وظيفي | CAD وB-Rep للتحقق من المسارين؛ اختبار sled لقياس ردود الفعل اليسرى واليمنى؛ مواءمة قوة المقعد مع FMVSS 207 وUNECE R17 |
| FR-02 | السماح بإزاحة طولية مضبوطة للمقعد عبر مسار ride-down `130` مع تبديد طاقة | V7 صفحات 3 و5 | أداء ديناميكي | علاقة force–stroke، hysteresis، temperature، cycling؛ بوابة V2 ثم V4 sled |
| FR-03 | فصل مسار دوران مسند الظهر ميكانيكياً عن مسار الإزاحة الطولية عبر `150L/150R` | V7 صفحات 3 و5 | استقلال وظيفي | اختبار carriage locked ثم translation permitted؛ قياس `theta(t)` وردود الوصلات؛ تحليل التداخل والـcoupling |
| FR-04 | انتقال حالات ميكانيكية يعتمد على حدث أو حالة حركة، ويشمل armed وride-down وrebound-control وsecure | V7 صفحة 3، V7 صفحة 5 | منطق حالة/سلامة وظيفية | اختبار تكرارية الانتقال، latency، عدم التحرير غير المقصود، fail-safe؛ لا يكفي إثباته بمخطط حالة فقط |
| FR-05 | التحكم في الحركة العكسية بعد الشوط الأساسي بواسطة `180` باستخدام damping أو friction أو locking | V7 صفحات 4 و5 | أداء/سلامة | قياس reverse travel وengagement point والسرعة والطاقة؛ اختبار مستقل ثم sled |
| FR-06 | بقاء مسار التقييد `190` قابلاً للعمل أثناء ride-down | V7 صفحة 3، V7 صفحة 5 | سلامة واجهة | تحليل مسار الحزام/التثبيت، اختبار anchorage، وحماية من التعارض أثناء الحركة؛ FMVSS 210 أو UNECE R14 بحسب سوق التطبيق |
| FR-07 | التحكم في حركة الحوض وتقليل forward migration عبر `140` مع توزيع الحمل على pan | V7 صفحات 4 و5 | تفاعل occupant/seat | sled وATD instrumentation؛ لا يثبت من الشكل الهندسي وحده |
| FR-08 | خرطوشة قابلة للاستبدال ذات علاقة force–stroke ومناطق مقاومة متعددة | V7 صفحة 4 | قابلية اختبار/صيانة | توصيف مستقل للخرطوشة، traceability، cycling، temperature، وتحقق أن التغيير لا يغير الهيكل الحامل |
| FR-09 | إمكانية قياس أو نمذجة ردود الفعل اليسرى واليمنى بصورة منفصلة | V7 صفحة 4 | قابلية القياس | قناتان مستقلتان، معايرة مشتركة وتزامن زمني؛ ISO 6487/SAE J211/1 كإطار قياس |
| FR-10 | تمثيل كل عنصر مرجعي كجزء حقيقي مع joint وtravel limit وmaterial definition وtest variable | V7 صفحة 9 | قابلية تنفيذ/توثيق | native parametric CAD، BOM، PMI/GD&T، material cards، evidence graph؛ الحالة الحالية BLOCKED |
| FR-11 | الحفاظ على stable reference numerals عبر CAD والرسومات والاختبارات وFEA | V7 صفحة 2 | traceability | تحقق آلي من mapping بين 100–210 والـBOM والاختبارات والـCAE |
| FR-12 | بوابات تحقق متدرجة من CAD إلى MBD والخرطوشة والقفل وsled وFE/HBM ثم full vehicle | V7 صفحة 8 | خطة تحقق | لا يجوز تجاوز بوابة إلى الادعاء النهائي؛ كل بوابة تحتاج artifact وhash وresult criteria |

## 4. المطابقة مع المعايير التنظيمية الحالية

### 4.1 FMVSS 207 — Seating Systems

يضع FMVSS 207 متطلبات لقوة المقعد، وتركيبه، واحتفاظه بوضع الضبط، وبعض المقاعد القابلة للطي. ويتطلب النص الحالي، ضمن شروط التطبيق والاستثناءات، اختبار قوى طولية أمامية وخلفية مرتبطة بكتلة المقعد، كما يتضمن عزماً عند الجزء العلوي من مسند الظهر. وعندما يكون حزام الأمان مثبتاً بالمقعد، تُطبق قوى الحزام بالتزامن مع قوى المقعد وفق FMVSS 210 [3].

| علاقة PDF | حالة المطابقة | المطلوب لإغلاق العلاقة |
|---|---|---|
| `100` و`110L/110R` و`120` كمسار حمل وتركيب للمقعد | مواءمة وظيفية قوية | تحديد فئة المركبة، اتجاه المقعد، كتلة المقعد، نقاط التثبيت، وحالات الضبط؛ تنفيذ اختبار FMVSS 207 المطبق وتوثيق القوة والتشوه والاحتفاظ بالضبط |
| `150L/150R` و`170` كتحكم مستقل لمسند الظهر والقفل | مواءمة جزئية | إثبات أن القفل لا يحرر أو يفشل تحت الأحمال المطلوبة، وأن آلية الطي/الضبط تحقق شروط القفل ذات الصلة |
| `130` ride-down كإزاحة مقصودة | لا توجد مطابقة تلقائية | يجب تحديد ما إذا كان المقعد يظل ضمن تعريف المقعد ومجال التطبيق، وكيف يحافظ على التثبيت أثناء الإزاحة؛ يلزم اختبار تنظيمي لا مجرد إظهار الشوط |
| `190` restraint interface | يحتاج FMVSS 210 متزامناً | لا يعالج FMVSS 207 وحده كل متطلبات التثبيت؛ يجب تحليل حالة belt-on-seat مع الأحمال المشتركة |

**حدود مهمة:** FMVSS 207 ليس معياراً لإثبات أن ride-down آمن أو أن claims الخاصة بـHIC/Nij محققة. كما أن إجراء NHTSA `TP-207-09` إرشادي للتنفيذ والتسجيل، وليس بديلاً عن نص اللائحة [4]. ويجب تحديد سنة تصنيع المركبة وفئة المركبة، لأن صياغة بعض الاستثناءات الحالية تغيرت في 2026 [5].

### 4.2 FMVSS 210 — Seat Belt Assembly Anchorages

يغطي FMVSS 210 نوع نقاط تثبيت الحزام، ومواقعها، وقوتها، وطريقة تحميلها، مع فروع مرتبطة بنوع الحزام وسنة تصنيع المركبة. ويتضمن النص الحالي فئات قوة مقدارها `5,000 lb` و`3,000 lb`، ومهلة تحميل مقدارها 10 ثوانٍ بعد الوصول إلى القوة، ومتطلبات هندسة الحزام ووضع الجسم الاختباري أو FAD [6].

| علاقة PDF | حالة المطابقة | المطلوب لإغلاق العلاقة |
|---|---|---|
| `190` occupant-restraint interface | مواءمة مباشرة من ناحية الوظيفة | released anchorage drawings، تحديد نوع الحزام، موضع التثبيت، seat/frame attachment، method المختار، وقوى الاختبار والتشوه |
| بقاء التقييد أثناء ride-down | غير مثبت | اختبار التثبيت مع مسار المقعد المتحرك، وإثبات عدم تعارض belt geometry أو زيادة slack أو تحميل غير مقصود |
| `210` sensor/trigger interface | ليست بديلاً عن anchorage evidence | فصل مستشعرات الحدث عن إثبات قوة التثبيت؛ يجب ألا يعتبر trigger signal دليلاً على صلاحية المرساة |

لا يجوز استخدام `18–22 kN` أو `25–35 kN` كبديل عن قوى FMVSS 210. تلك القيم أهداف مصدرية للـabsorber، بينما FMVSS 210 يحدد أحمالاً تنظيمية مختلفة بحسب نوع التثبيت والحالة القانونية.

### 4.3 FMVSS 202a — Head Restraints

يحدد FMVSS 202a متطلبات توفير مسند الرأس، وارتفاعه وعرضه والفجوات والاحتفاظ بالضبط والقوة وامتصاص الطاقة. كما يوفر مساراً ديناميكياً يتضمن حد دوران خلفي بين الرأس والجذع مقداره 12 درجة و`HIC15` مقداره 500، مع شروط الاختبار الأخرى [7].

| علاقة PDF | حالة المطابقة | المطلوب لإغلاق العلاقة |
|---|---|---|
| `160` seatback frame و`200` torso/head guidance | مواءمة وظيفية جزئية | تحديد ما إذا كان `200` head restraint بالمعنى التنظيمي أم guidance geometry مختلفة؛ قياس H-point، torso reference line، backset، height، width، gap، strength، energy absorption |
| التحكم في trajectory | أوسع من نطاق FMVSS 202a | لا يجوز اعتبار كل توجيه للرأس أو الجذع مطابقاً لمعيار head restraint؛ يجب فصل المتطلبات التنظيمية عن فرضيات trajectory control |
| seatback rotation path | تداخل محتمل | اختبار head-restraint في أوضاع الضبط والأسوأ، مع تحليل تأثير دوران المسند على المرجع الهندسي والـbackset |

المعيار لا يثبت وحده أن pelvic control أو rebound control آمنان. كما أن مسار 12 درجة/HIC15 ليس تصريحاً بأن النظام يحقق حماية إصابات شاملة؛ يجب استيفاء باقي شروط المسار المختار.

## 5. المطابقة مع UNECE Regulations

| المعيار | الوظيفة التي يغطيها | الصلة بـV7 | النتيجة الحالية |
|---|---|---|---|
| UN R14 | تثبيتات أحزمة الأمان ومواقعها وقوتها، بما فيها التثبيت على بنية المقعد في الحالات ذات الصلة | `190` ومسار الحمل بين seat وvehicle | مواءمة مباشرة، لكن لا توجد رسومات أو test evidence مفرج عنها |
| UN R17 | قوة المقاعد وتثبيتات المقاعد، الضبط والقفل، وبعض متطلبات مساند الرأس والحماية الخلفية | `100` و`120` و`160` و`170` | مواءمة وظيفية جزئية؛ يلزم تحديد فئة المركبة واتجاه المقعد وسلسلة التعديل المطبقة |
| UN R25 | اعتماد جهاز مسند الرأس عندما لا يغطيه مسار R17 | `200` إذا كان head-restraint device مستقلاً | مجال محتمل فقط؛ لا يجوز تطبيق R25 وR17 معاً بلا تحديد قاعدة عدم الازدواج |

تنص مواد UNECE على فصل الوظائف: R14 للحزام، R17 للمقعد وتثبيته وبعض مساند الرأس، وR25 لجهاز مسند الرأس عندما لا يكون مسار R17 هو المسار الحاكم. لذلك فإن claim المعماري الواحد في PDF يحتاج **مصفوفة اعتماد متعددة** بدلاً من الادعاء بأن معياراً واحداً يغطي النظام كله [8] [9] [10].

## 6. القياس والمعايرة: ISO 6487 وSAE J211/1

تدعم ISO 6487:2015، مع تعديل 2017، وSAE J211/1_202208 بناء قناة القياس كاملة من الحساس والتوصيلات إلى التحليل. وتشمل الموضوعات ذات الصلة: `CFC`، استجابة التردد والسعة، phase delay، time base، relative time delay، المعايرة، عدم الخطية، الحساسية، وعدم اليقين [11] [12].

هذا يطابق احتياجات PDF المتعلقة بـ`210`، وبالقياس المستقل لردود الفعل اليسرى واليمنى، وبوابات sled وFE/HBM. لكنه لا يحدد وحده قيمة نجاح أو فشل للمقعد أو الإصابة. يجب أن يأتي معيار القبول من FMVSS أو UNECE أو بروتوكول اختبار أو عقد اختبار محدد.

ويجب الانتباه إلى أن ISO 6487 تحت المراجعة، وأن النصوص الكاملة محمية بحقوق النشر. لذلك لا ينبغي تثبيت أرقام CFC أو حدود زمنية من نسخ قديمة أو previews دون التحقق من الإصدار الحالي المطلوب تعاقدياً.

## 7. المتطلبات التي لا تغطيها هذه المعايير مباشرة

بعض المتطلبات الوظيفية في PDF لا تملك مقابلاً تنظيمياً مباشراً في المعايير المذكورة. وتشمل ذلك منطق الحالات `armed → ride-down → rebound-control → secure`، وفصل مسار دوران الظهر عن ride-down كعلاقة ابتكارية، وخوارزمية الانتقال المبنية على التسارع أو السفر أو الدوران، وخرطوشة force–stroke ذات منطقتي مقاومة، وقابلية إعادة استخدام المنطق في تغليف أمامي وخلفي.

هذه المتطلبات تحتاج مواصفات هندسية داخلية قابلة للقياس. ويجب أن تتضمن كل مواصفة تعريف الحالة، الشرط الابتدائي، trigger، latency، القوة أو الحركة المتوقعة، حالة الفشل، معيار القبول، وعدم التحرير غير المقصود. لا يجوز اعتبارها مغلقة بالاستناد إلى FMVSS 207 أو R17 وحدهما.

## 8. مصفوفة الأدلة المطلوبة للإطلاق الهندسي

| طبقة الدليل | الحد الأدنى المطلوب | حالة V7 الحالية |
|---|---|---|
| CAD/B-Rep | native parametric assembly، joints، travel limits، interference، envelope، reference numerals | **MISSING/BLOCKED** |
| BOM والمواد | released BOM، material grade، certificates، density، lot traceability | **MISSING/BLOCKED** |
| PMI/GD&T | datums، dimensions، tolerances، hole patterns، clearances، inspection method | **MISSING/BLOCKED** |
| interfaces | vehicle hardpoints، seat anchorage، belt path، adjustment and locking interfaces | **MISSING/BLOCKED** |
| absorber | force–stroke، hysteresis، temperature، cycling، serial identity | **MISSING/BLOCKED** |
| lock/rebound | state timing، loads، unintended-release، fail-safe، reverse-travel evidence | **MISSING/BLOCKED** |
| instrumentation | sensor ID، calibration، CFC/filter، timebase، synchronization، uncertainty، raw-data hash | **PARTIAL CONTRACT ONLY** |
| sled/ATD | calibrated channels، high-speed kinematics، seat loads، posture، pre/post inspection | **NOT AVAILABLE** |
| full vehicle | vehicle-level regulatory/consumer evidence | **NOT AVAILABLE** |

## 9. القرار والحدود

الوثائق المرفقة صالحة لتحديد **functional engineering requirements** ولإنشاء verification plan، لكنها لا تكفي لإعلان compliance أو manufacturability أو crash-safety performance. المطابقة الحالية هي مواءمة تحليلية جزئية مع المعايير، وليست شهادة.

يجب إبقاء الادعاءات التالية مفتوحة: أن `180 mm` آمن أو محقق، أن `18–22 kN` يصف absorber فعلياً، أن `25.6 g` هو pulse معتمد للنظام، أن HIC/Nij سيتحسنان، وأن المعمارية تحقق FMVSS أو UNECE. إغلاق أي منها يتطلب مصدر authoritative أو اختباراً/تحليلاً موثقاً مع provenance وhash وconfiguration mapping.

بناءً على ذلك، يظل قرار المشروع: **STOP RELEASE**. الخطوة التالية المقبولة هندسياً هي استرداد V7 native package أو اعتماد بديل رسمي من صاحب الصلاحية، ثم بناء compliance matrix مرتبطة بفئة المركبة وسوقها وسنة التصنيع وبالإصدار التنظيمي المحدد.

## References

[1]: file:///home/ubuntu/upload/T_OCS_V7_Patent_Draft_Invention_Disclosure_.pdf "T-OCS V7.0 Patent Preparation — V7 Enhancement Addendum"

[2]: file:///home/ubuntu/upload/نماذجV5-V7.pdf "T-OCS V5–V7 Engineering Forms and Visual Boards"

[3]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.207 "49 CFR 571.207 — Standard No. 207; Seating systems"

[4]: https://www.nhtsa.gov/sites/nhtsa.gov/files/2023-06/TP-207-09.pdf "NHTSA OVSC TP-207-09 — FMVSS 207 Seating Systems"

[5]: https://www.federalregister.gov/documents/2026/06/03/2026-11078/advanced-improvements-in-school-bus-occupant-protection "Federal Register 91 FR 33098 — Advanced Improvements in School Bus Occupant Protection"

[6]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.210 "49 CFR 571.210 — Standard No. 210; Seat belt assembly anchorages"

[7]: https://www.ecfr.gov/current/title-49/subtitle-B/chapter-V/part-571/subpart-B/section-571.202a "49 CFR 571.202a — Standard No. 202a; Head restraints"

[8]: https://unece.org/sites/default/files/2025-05/R014r6e%20%282%29.pdf "UNECE UN Regulation No. 14 — Safety-belt anchorages"

[9]: https://unece.org/sites/default/files/2022-01/R017r6e.pdf "UNECE UN Regulation No. 17 — Seats, their anchorages and head restraints"

[10]: https://unece.org/sites/default/files/2025-08/R025r1a4e.pdf "UNECE UN Regulation No. 25 — Head restraints"

[11]: https://www.iso.org/standard/64041.html "ISO 6487:2015 — Measurement techniques in impact tests"

[12]: https://saemobilus.sae.org/standards/j2111_202208-instrumentation-impact-test-part-1-electronic-instrumentation "SAE J211/1_202208 — Instrumentation for Impact Test, Part 1: Electronic Instrumentation"
