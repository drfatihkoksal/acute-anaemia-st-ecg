## 11. PROTOKOL v3 — KİLİTLİ (2026-10-02): HEEDB'de keşif, MIMIC'te dış doğrulama

**Durum:** §1–§10 tarihsel kayıttır. Psödo-iskemi/düzeltici-baş çerçevesi bu protokolle **değiştirilmiştir**. Bu bölüm, MIMIC (BIDMC) verisi analiz amacıyla açılmadan önce kilitlenmiştir. Kilit sonrası her değişiklik §11.9 değişiklik günlüğüne tarih ve gerekçeyle yazılır ve makalede raporlanır.

### 11.1 Gerekçe (önceki bulgular + fizibilite)
- MIMIC: çift-EKG ΔHb r≈0.25 (tek EKG 0.15, yanlış hasta çifti 0.12). Kreatinin kontrolü de küçük bir kazanç gösterdi (0.08). Psödo/gerçek iskemi ayrımı başarısız (AUC≈0.5, n_gerçek=40).
- HEEDB-MGH fizibilitesi (`heedb/01`, `heedb/02`): 898.761 çift / 172.436 hasta. 12SL iskemik ifade oranı Hb düşüşünden sonra artmıyor (%15.3; stabil Hb'de %17.2).
- Literatür (`heedb/literature_H2_brody.md`):
  - Akut izovolemik anemide ST değişikliği nadir ve kalp hızıyla ilişkili (Weiskopf 1998; Leung 2000).
  - Hct düşüşü radyal (erken) QRS'i büyütüp tanjansiyel (geç) QRS'i küçültüyor (Brody etkisi; Rosenthal 1971, Nelson 1972, Hodgkin 1977, Oreto 1992). Ancak tutarsız (Ishikawa 1976).
  - Volüm ve ödem etkileri bütün QRS ile P'yi aynı yönde ölçekliyor (Madias 2005, Drighil 2006, Kataoka 2011).

### 11.2 Araştırma sorusu
Hasta kendi kontrolü olacak şekilde bakıldığında, akut Hb düşüşü EKG'yi iskemik bir ST-T değişikliğiyle mi, yoksa kalp hızından bağımsız, kan iletkenliğine bağlı (Brody) bir QRS amplitüd örüntüsüyle mi değiştirir?

### 11.3 Hipotezler (yön ve karar kuralı önceden belirlenmiştir)
Maruziyet: ΔHb = Hb_baz − Hb_son (g/dL; pozitif = düşüş). Tüm etkiler, 1 g/dL düşüş başına β olarak raporlanır.

- **H1, ST eşdeğerliği (birincil):**
  - Sonuç: ΔST60_min = 12 derivasyondaki en düşük ST60 (QRSoff+60 ms, PR bazal çizgisine göre) değerinin değişimi.
  - Hipotez: Kovaryat düzeltmesinden sonra |β| < 0.01 mV/(g/dL).
  - Gerekçe: 3 g/dL düşüş için 0.03 mV eder; bu, 0.05 mV'lik klinik ST depresyonu eşiğinin altında kalır.
  - Test: TOST. 90% CI tamamen ±0.01 içindeyse eşdeğerlik kabul edilir. Sonuç, MGH ve BIDMC'nin **ikisinde de** sağlanmalı.
  - İkincil sonuçlar: Derivasyon başına ΔST60, ΔJ, ΔT amplitüdü. Ayrıca kalp hızı düzeltmesi **olmadan** aynı model (taşikardi aracılığının büyüklüğünü göstermek için).
- **H2a, Brody parmak izi (birincil):**
  - Kors dönüşümüyle VKG. Medyan atımda QRSon–QRSoff aralığı normalize edilir:
    - erken pencere = %0–45, en yüksek uzaysal büyüklük M_erken;
    - geç pencere = %55–100, en yüksek uzaysal büyüklük M_geç.
  - Öngörü: Hb düşüşüyle β(ΔM_erken) > 0 ve β(ΔM_geç) < 0.
  - Karar: İki β'nın 95% CI'ı sıfırı dışlıyor ve işaretler zıt ise, MGH'de "destekleniyor" denir. BIDMC'de işaretler aynıysa ve CI'lar sıfırı dışlıyorsa "replike" denir.
  - İkincil (Oreto replikasyonu):
    - ΣR (12 derivasyon);
    - S_V6 = V6'da QRS penceresindeki en düşük değer (bazal çizgiye göre; S yoksa 0).
    - Öngörü: Hb düşüşüyle ΣR ↑ ve |S_V6| ↓.
- **H2b, rakip volüm-iletken hipotezi:**
  - Öngörü: ΣQRS (tepe-tepe, 12 derivasyon) ile ΣP aynı yönde değişir. ΔM_erken ile ΔM_geç de aynı işaretli olur.
  - Bu ortak bileşen, ΔHb'den çok volüm göstergeleriyle ilişkilidir: ΔNa, ΔBUN, Δalbümin (hangisi mevcutsa).
  - Karar: H2a reddedilir ve ortak ölçekleme anlamlıysa, sonuç "volüm baskın" olarak raporlanır.
- **H3, bilginin yeri (ikincil, derin öğrenme):**
  - Mevcut çift-füzyon modeli üç girdi koşuluyla eğitilir: yalnız QRS penceresi, yalnız ST-T penceresi, tam sinyal.
  - Öngörü: r(yalnız QRS) ≥ 0.8 × r(tam sinyal) ve r(yalnız ST-T) < 0.5 × r(tam sinyal).
  - Eğitim yalnızca MGH'de, değerlendirme BIDMC'de yapılır.
- **H4, dış geçerlilik:** H1–H3'ün tüm tahminleri, model veya eşik yeniden ayarlanmadan BIDMC'de raporlanır. Merkezler arası fark için heterojenlik testi yapılır.

### 11.4 Kohortlar
- **Keşif:** HEEDB I0001 (MGH), `heedb/out/pairs_mgh.parquet`.
  - Her EKG ±12 saat içindeki en yakın Hb'ye eşlenir (LOINC 718-7, 30313-1, 30350-3, 20509-6; durum Final/Corrected; 2–25 g/dL).
  - Ardışık iki EKG, 0 < aralık ≤ 7 gün ise ve farklı Hb ölçümlerine bağlıysa çift oluşturur.
  - Not: HEEDB EKG metadata'sında yatış bağlantısı yoktur. MIMIC'teki "aynı yatış" şartı yerine burada "Hb ±12 saat" şartı vekil olarak kullanılır.
- **Dış doğrulama:** MIMIC-IV-ECG + MIMIC-IV (BIDMC), `scripts/30` tanımıyla `manifests_aligned`.
  - Doğrulama, aynı yatış şartı **dahil** orijinal tanımla yapılır.
  - Duyarlılık analizi: yatış şartı kaldırılarak HEEDB tanımıyla tekrarlanır.
- **Hasta çakışması:** HEEDB'de BIDMC verisi de vardır. MGH ile BIDMC arasında hasta bağlantısı kurulamaz ve kurulmayacaktır; bu yüzden olası ortak hastalar sınırlılık olarak yazılır.
- **Hariç tutma (iki merkezde aynı, sinyalden türetilmiş):**
  - çiftin herhangi bir EKG'sinde ölçüm başarısızlığı (ok=False; tutulan atım < 6);
  - ΔQRS süresi > 20 ms (yeni dal bloğu veya iletim değişikliği vekili);
  - RR düzensizliği ritim değişikliği gösteriyorsa (çiftin birinde RR değişim katsayısı > 0.15, diğerinde değilse);
  - pace ritmi: HEEDB'de 12SL kodları, MIMIC'te makine raporu, ayrıca sinyalde pace spike'ı.
- **Birincil analiz birimi:** Hasta başına bir çift; deterministik hash ile seçilir. Duyarlılık analizinde tüm çiftler kullanılır, standart hatalar hastaya göre kümelenir.

### 11.5 Ölçüm hattı (iki merkezde aynı, dondurulmuş) — v3.1 (2026-10-02)
- **Modül:** Wavedet kullanılmıyor (bkz. §11.9). Ölçüm `heedb/ecgmeasure_gpu.py` (sha1 f47a1270) ile GPU'da (PyTorch) yapılır; toplu çalıştırıcı `heedb/04_measure.py` (sha1 de33fc0b).
- **Adımlar:**
  1. 0.5–100 Hz sıfır-fazlı FFT filtresi (filtfilt eşdeğeri Butterworth) ve 60 Hz notch.
  2. 12 derivasyonun türev enerjisinden QRS tespiti (250 ms refrakter).
  3. Ektopik atım kuralı: RR < 0.8 × medyan RR ise o atım ve sonraki atılır.
  4. Atımlar ±20 ms çapraz korelasyonla hizalanır; şablon korelasyonu ≥ 0.90 olanlar tutulur.
  5. Medyan atım: −300 ms ile +500 ms arası.
  6. QRS başlangıç/bitişi 12 derivasyon uzaysal hızından bulunur. Eşik gürültüye uyarlanır: QRS'in öncesindeki ve sonrasındaki en düz 20 ms'lik pencerenin gürültü seviyesi + a × (tepe hız − gürültü).
  7. Bazal çizgi: QRSon'dan −40 ile −12 ms önceki aralığın medyanı.
  8. Kors VKG; özellikler §11.3'teki gibi.
- **Eşik ayarı (LUDB, `heedb/03_validate_ludb.py`):** Ayar tek numaralı kayıtlarda, test çift numaralı kayıtlarda yapıldı.
  - Seçilen değerler: a_on = a_off = 0.05; sabit sapma düzeltmesi başlangıç için −14 ms, bitiş için +14 ms.
  - Test sonucu (n=98):
    - başlangıç: ortalama −0.3 ms, SD 8.8 ms, MAE 6.1 ms, %89'unda hata ≤10 ms;
    - bitiş: ortalama −1.7 ms, SD 17.9 ms (uç değerler), MAE 8.5 ms, %79'unda ≤10 ms.
  - Seçilen a, grid'in üst sınırına denk geldi.
- **QRS sınır tanımı:** Sınırlar LUDB'nin "küresel" tanımını izler (12 derivasyonda en erken başlangıç → en geç bitiş). LUDB'de medyan QRS süresi 121 ms. MIMIC makine QRS süresinden medyan 22 ms daha uzun, aralarındaki korelasyon r=0.83. ST60 bu J noktasına göre ölçülür; makaleye bu tanımla yazılır.
- **Hariç tutma:** Tutulan atım sayısı < 6 ise kayıt başarısız sayılır (§11.4'teki 8'in yerine; bkz. §11.9).
  - Başarı oranı: MGH 1.365.629 EKG'de %96.2; MIMIC 137.973 EKG'de %94.8 (311'i NaN sinyal).
- **Örnekleme:** 250 Hz kayıtlar `resample_poly` ile 500 Hz'e çıkarılır. Duyarlılık analizinde yalnız 500 Hz kayıtlar kullanılır.
- **Kör ölçüm:** Ölçüm hattı ΔHb'yi görmez. MIMIC özellikleri hesaplandı ama ΔHb ile hiçbir analiz yapılmadı (§11.8 geçerli).

### 11.6 İstatistik modeli
- Birincil model: Δözellik ~ ΔHb + ΔRR + ΔK + ΔCa + aralık_saat + yaş + cinsiyet + baz_özellik.
- Elektrolitler: EKG'ye ±12 saat içinde en yakın değer. Eksikse çoklu atama yapılır; duyarlılık analizi tam veri ile.
- ΔHb'nin doğrusal olmayan etkisi kübik spline ile ikincil analizde incelenir.
- Çoklu test: Birincil sonuçlar H1, H2a(M_erken), H2a(M_geç) için Holm düzeltmesi uygulanır. İkincil sonuçlar keşifsel olarak etiketlenir.

### 11.7 Negatif kontroller ve alt gruplar
- **Stabil Hb çiftleri (|ΔHb| < 0.5):** Ölçüm ve elektrot gürültüsünün tabanını verir. Etkiler bu tabana göre yorumlanır.
- **Yanlış hasta çiftleri:** Etki ortadan kalkmalı.
- **Alternatif maruziyet ΔKreatinin:** Brody parmak izini üretmemeli.
- **Alt gruplar (önceden tanımlı):**
  - kanama vekili: GI kanama veya hemoraji ICD kodu ±3 gün, ya da transfüzyon kaydı;
  - hemodilüsyon/postop vekili: kanama kodu yok ve cerrahi kodu var, ya da Δalbümin ve ΔNa düşüşü;
  - Hb ≤ 8 olan son EKG'ler.
  - Öngörü: Brody parmak izi iki alt grupta da görülür. Volüm bileşeni alt gruplar arasında farklı yönde olabilir.

### 11.8 Yayın ve süreç kuralları
- MIMIC sonuçlarına, MGH analiz kodu dondurulup hash'i kaydedilmeden bakılmaz.
- H1 için her iki sonuç ("eşdeğer" ya da "eşdeğer değil") raporlanır. H2a reddedilirse negatif sonuç olarak raporlanır.
- **Veri kullanım sözleşmeleri:** BDSP (HEEDB) ve PhysioNet (MIMIC) şartlarının, aynı makinede iki veri setiyle analiz yapmayı ve model ağırlıklarını yayımlamayı kapsadığı, analizden önce teyit edilecek. Hastalar arası bağlantı veya eşleştirme yapılmaz.

### 11.9 Değişiklik günlüğü
- 2026-10-02: v3 kilitlendi.
- 2026-10-02 (v3.1, MIMIC'e bakılmadan önce): Kullanıcı kararıyla ölçüm hattı wavedet'ten GPU tabanlı medyan atım + VKG modülüne geçirildi (§11.5). Gerekçeler: bağımlılıksız Python, ~8.000 EKG/sn (wavedet ~15), protokoldeki VKG ölçümleriyle doğrudan uyum. Asgari atım sayısı 8'den 6'ya indirildi: 10 sn kayıtta HR 50'de ~8 atım olur, ektopik kuralı 2 atım siler. Hipotezler, karar kuralları ve kohort tanımları değişmedi.
- 2026-10-02: Kors matrisi BRAVEHEART'ın açık kaynak `korsmatrix.m` dosyasıyla (commit 9b2e03b) karşılaştırıldı: aynı, en büyük fark 2e-8. Not: Kors 1990'ın orijinal tablosuyla doğrudan değil, ikincil bir kaynakla karşılaştırıldı.
- 2026-10-02: MGH keşif analizi çalıştırıldı (`heedb/05_build_analysis_mgh.py` sha1 1e05bd4f, `heedb/06_analysis.py` sha1 788c5b4d). **Analiz kodu bu hash'lerle donduruldu.** MIMIC analizi `06_analysis.py` değiştirilmeden çalıştırılacak; değişmesi gereken tek parça MIMIC'e özgü veri inşa betiği.
  - Uygulama notları:
    - Pace bayrağı MGH'de 12SL kodlarından alındı.
    - Elektrolit eksikliği %5.2; çoklu atama ile (m=5) dolduruldu.
    - Analiz akışı: 898.761 çift → 663.991 çift → hasta başına bir çiftle 154.187 hasta.
  - **Birincil sonuçlar (MGH):**
    - **H1 eşdeğerliği sağlandı:** β(ΔST60_min) = +0.0005 mV/(g/dL), 90% CI [0.0004, 0.0007]; TOST p<1e-15.
    - **H2a reddedildi:** M_erken β = −0.0052 [−0.0064, −0.0040], M_geç β = −0.0109 [−0.0121, −0.0096]. İki bileşen aynı yönde küçülüyor.
  - **İkincil bulgular:**
    - ΣQRS ve ΣP birlikte azalıyor; bu azalma Δalbümin ve ΔNa ile güçlü ilişkili (H2b ile uyumlu).
    - Kreatinin maruziyeti de aynı örüntüyü üretiyor (özgül değil).
    - Hb düşüşü kalp hızının **azalmasıyla** birlikte görülüyor (β −0.80 vuru/dk per g/dL). Bu, maruziyetin büyük ölçüde hemodilüsyon/volüm yerine koyma tablosunu yansıttığına işaret ediyor.
    - Kanama ICD alt grubunda (n=3.262) yön Brody ile uyumlu (M_erken +0.005, p=0.056; M_geç −0.003, anlamsız), ama güç yetersiz. Sonuç keşifsel olarak raporlanacak.
- 2026-10-02 (MIMIC analizi çalıştırılmadan önce): MIMIC inşa betiği `heedb/05b_build_analysis_mimic.py` (sha1 ffff2916) yazıldı. `06_analysis.py` değişmedi (sha1 788c5b4d).
  - MIMIC'e özgü zorunlu farklar:
    - kanama ICD kodu "aynı yatış" düzeyinde (MIMIC'te tanılar tarihsiz). Oranı %24.6, MGH'de ±3 gün tanımıyla %1.9; bu yüzden kanama alt grubu iki merkez arasında karşılaştırılabilir değil;
    - pace bayrağı makine raporu metninden (regex) alındı;
    - lab itemid'leri MGH'deki LOINC setinin karşılıkları; albümin kapsamı %21.
- 2026-10-02: **MIMIC dış doğrulaması** (donmuş kod, sha1 788c5b4d). Akış: 87.867 çift → 58.506 çift → hasta başına bir çiftle 24.741 hasta.
  - **H1 replike oldu:** β = +0.0004 mV/(g/dL), 90% CI [−0.0001, 0.0008], TOST p≈0. Kalp hızı düzeltmesi olmadan da sıfır.
  - **H2a yine reddedildi:** M_erken −0.0169 [−0.0197, −0.0141], M_geç −0.0192 [−0.0220, −0.0165]. İki bileşen aynı yönde küçülüyor. Etki MGH'dekinin 2–3 katı; heterojenlik anlamlı (M_erken z≈−7.5).
  - **Volüm bileşeni:** ΣQRS ve ΣP'nin birlikte azalması ve albümin ilişkisi replike oldu. ΔNa ilişkisinin işareti ters döndü. Volüm düzeltmesi MIMIC'te erken QRS etkisini kaldırmadı (albümin kapsamı %21).
  - **Diğer:**
    - Kreatinin maruziyeti MIMIC'te Hb'ye göre zayıf; yani bu merkezde QRS küçülmesi daha Hb'ye özgü.
    - Yanlış hasta kontrolleri sıfır.
    - Kanama alt grubunda (yatış düzeyinde tanım) iki bileşen de küçülüyor. MGH'deki Brody yönündeki keşifsel ipucu replike olmadı; ancak alt grup tanımları karşılaştırılabilir değil.
- 2026-10-02 **KEŞİFSEL (önceden kayıtlı değil)**, `heedb/07_exploratory.py` (sha1 ae5d3638). Fonksiyonlar donmuş 06'dan alındı; tüm tahminler tam veri (complete-case).
  - **Heterojenlik:** Tam veride erken QRS β MGH'de −0.0054, MIMIC'te −0.0123; geç QRS −0.0111 ve −0.0163. Bu analizlerle heterojenlik açıklanamadı:
    - MGH'yi aynı encounter'a ait çiftlerle sınırlamak (%77) MIMIC'e yaklaştırmadı (erken −0.0044).
    - MGH'yi MIMIC'in kovaryat dağılımına ağırlıklandırmak da yaklaştırmadı (erken −0.0028).
    - Tanımlayıcı özellikler iki kohortta benzer.
    - Kalan adaylar: cihaz/filtre farkları, Hb ölçüm yöntemi, ölçülmemiş sıvı uygulamaları.
  - **Kesin kanama (D62 / ICD-9 285.1):** MGH'de yalnız n=142 (D62 nadiren kodlanmış), MIMIC'te n=3.501.
    - ST: iki merkezin birleşik tahmini −0.0005 [−0.0017, 0.0008]; heterojenlik p=0.66. Eşdeğerlik marjının içinde, yani **H1 kesin kanamada da geçerli**.
    - QRS: birleşik erken QRS −0.0067 [−0.0140, 0.0007], geç QRS −0.0069 [−0.0129, −0.0009]. Brody örüntüsü yok.
  - **Majör kanama kodları:** MGH n=2.818 (erken QRS +0.0047 ve ΣQRS ≈0, ikisi de anlamsız), MIMIC n=4.724 (erken −0.0080, geç −0.0093, ΣQRS −0.12).
- 2026-10-02 **KEŞİFSEL**, `heedb/08_tail_st.py` (sha1 88b63faf). Soru (kullanıcı itirazı): "Ortalama sıfır olsa da düşük Hb'de bir alt grupta ST depresyonu ve anemiye bağlı MI gelişiyor mu?"
  - Sonuç değişkeni: yeni ST depresyonu, yani aVR dışında en az 2 derivasyonda post ST60 ≤ −0.05 mV ve ΔST60 ≤ −0.05 mV.
  - **Arka plan oranı** (Hb stabil, nadir >10): MGH %5.4, MIMIC %6.9. Bu, ölçüm gürültüsü ve diğer nedenlerin payını gösteriyor.
  - **Düşük nadir Hb ile oran artıyor:**
    - Hb ≥1.5 düşüp ≤7'ye inenlerde: MGH %9.6 (n=904), MIMIC %13.5 (n=148). Düşüp >10'da kalanlarda: %5.4 ve %7.5.
    - Ancak Hb'si stabil olup ≤7 olanlarda da oran yüksek (MGH %7.3, MIMIC %8.8). Yani artışın bir kısmı düşüşten değil, düşük Hb düzeyinden ve hastalık ağırlığından geliyor.
  - **Lojistik regresyon:**
    - Hb ≤8 OR: MGH 1.18 [1.10, 1.26], MIMIC 1.24 [1.08, 1.43].
    - ΔHb tek başına OR ≈1.
    - ΔHb × nadir ≤8 etkileşimi: MGH'de 1.19 [1.13, 1.25]; MIMIC'te 1.06, anlamsız.
    - Kalp hızı değişimi modele eklenince ΔHb etkisi tamamen kayboluyor. Bu, taşikardi üzerinden bir aracılıkla uyumlu.
  - **Troponin** (ölçülme oranı %16–22, seçici ölçüm):
    - Troponin yükselmesi ΔHb ile artıyor: MGH %30 (stabil) → %39 (≥2.5 düşüş); MIMIC %66 → %79.
    - Troponin yükselenlerde yeni ST depresyonu daha sık: MGH %9.6 / %6.3; MIMIC %13.4 / %10.2.
  - **Yorum:** H1 (ortalama etki yok) korunuyor, ama kapsamı daraltılmalı. Ağır anemide (nadir ≤7–8) yeni ST depresyonu ve miyokard hasarı belirgin artıyor. Etki taşikardiyle aracılanıyor ve küçük bir azınlıkta yoğunlaşıyor. Makalede bu "kuyruk bulgusu" ikincil bulgu olarak ve H1 ile birlikte raporlanmalı.
- 2026-10-02 **Yeni soru için fizibilite** ("EKG ile anemi intoleransı"), `heedb/09_feasibility_tolerance.py` (sha1 ef634b15). Bu yeni bir çalışmanın ön sayımıdır, v3 protokolünün parçası değildir.
  - İndeks olay: Hb ≥1 g/dL düşüp ≤8 g/dL'ye inen ilk analiz edilebilir çift.
  - **MGH:** 18.288 hasta (Hb ≤7: 5.071). Yeni ST depresyonu %7.9 (n=1.443).
    - 30 günlük mortalite %16.9 (3.075 ölüm): ST depresyonu olanlarda %25.0, olmayanlarda %16.2.
    - 1 yıllık mortalite %33.0: %40.9'a karşı %32.4.
    - Troponin ölçülme oranı %29.
  - **MIMIC:** 1.928 hasta. Yeni ST depresyonu %11.6.
    - 30 günlük mortalite %18.2 (351 ölüm): %28.3'e karşı %16.9.
    - 1 yıllık mortalite %38.2.
  - Bu oranlar ham (düzeltilmemiş) değerlerdir; hastalık ağırlığına göre düzeltme gerekiyor.

---

