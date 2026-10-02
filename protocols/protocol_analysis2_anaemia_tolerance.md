# Protokol: Akut Anemide EKG ile Anemi İntoleransı (prognostik çalışma)

**Durum:** KİLİTLİ v1.0 (2026-10-02, kullanıcı onayı). Bölüm 9'daki üç karar önerildiği gibi kabul edildi. Kod MGH sonuç analizinden önce dondurulacak ve hash'leri §10'a yazılacak.
**Bağlam:** concept.md §11 (protokol v3) ve §11.9 günlüğü. Ön sayım: `heedb/09_feasibility_tolerance.py`. Ölçüm hattı: `heedb/ecgmeasure_gpu.py` (§11.5 v3.1, değişmeden kullanılır).

## 1. Gerekçe
- **Bizim bulgularımız:**
  - Hb düşüşü ortalamada ST'yi kaydırmıyor (v3 H1).
  - Nadir Hb ≤7–8 olduğunda yeni ST depresyonu ~%5–7'den ~%10–13'e çıkıyor ve bu artış taşikardiyle aracılanıyor (`08_tail_st.py`).
- **Literatür:**
  - Transfüzyon kararı Hb eşiğine dayanıyor. MINT çalışmasında kısıtlayıcı strateji kardiyak ölümü artırdı (meta-analiz RR 1.47).
  - Fizyolojik tetikleyiciler (ScvO2, O2ER, laktat) yoğun bakıma ya da invaziv ölçüme bağlı.
  - 12 derivasyonlu seri EKG ile anemi intoleransını sonlanımlara karşı doğrulayan büyük bir çalışma bulunamadı. En yakınlar: travmada kalp hızı değişkenliği ile makine öğrenmesi (Gopalakrishnan 2025) ve küçük bir kraniyotomi çalışması.
- **Ön sayım (ham oranlar):** MGH'de 18.288 hasta, 30 günlük mortalite yeni ST depresyonu olanlarda %25.0, olmayanlarda %16.2. MIMIC'te 1.928 hasta, %28.3'e karşı %16.9.

## 2. Araştırma sorusu
Hb'si akut olarak ≤8 g/dL'ye düşen erişkinlerde, önceki EKG'ye göre gelişen yeni ST depresyonu, Hb düzeyi ve hastalık ağırlığından bağımsız olarak 30 günlük mortaliteyi öngörür mü? Bu prognostik değer, anemisi olmayan hastalardaki genel ST depresyonu riskinden daha mı büyüktür?

## 3. Kohortlar
- **İndeks olay:** Hasta başına ilk analiz edilebilir EKG çifti; çift tanımı v3 §11.4 ile aynı. Şartlar: son EKG'deki Hb ≤ 8.0 g/dL ve ΔHb (baz − son) ≥ 1.0 g/dL. İndeks zamanı = son EKG.
- **Hariç tutma:** v3 `06_analysis.prep` kuralları (iki EKG de ölçülmüş, |ΔQRS| ≤20 ms, ritim değişikliği yok, pace yok) ve yaş <18.
- **Takip:**
  - MGH: Massachusetts ölüm kaydı (2023-04-30'a kadar tam). 30 günlük takibi tamamlanmamış indeksler çıkarılır.
  - MIMIC: patients.dod. Son taburculuktan 1 yıl sonrası sansürlü olduğu için 30 günlük takip tam kabul edilir; 1 yıllık analizde bu sansür belirtilir.
- **Kontrast kohortu (özgüllük testi):** Hb'si stabil (|ΔHb| < 0.5) ve > 10 g/dL olan çiftler; aynı kurallarla hasta başına ilk çift. Anemik kohortla aynı yıl ve yaş dağılımına gelecek şekilde ağırlıklandırılır.
- **Keşif:** MGH. **Dış doğrulama:** MIMIC (BIDMC).

## 4. Maruziyet
- **Birincil:** Yeni ST depresyonu. Tanım: aVR dışında en az 2 derivasyonda hem son EKG'de ST60 ≤ −0.05 mV, hem de ΔST60 ≤ −0.05 mV.
  - Not: ST60, LUDB "küresel" J noktasına göre ölçülür (v3 §11.5); bu J noktası cihazın J noktasından ~10–20 ms daha geçtir.
- **İkincil:**
  - Aynı tanımın 0.10 mV (1 mm) eşikli hali.
  - Sürekli "yeni ST depresyon yükü" = Σ max(0, −ΔST60) (0.05 mV kriterini karşılayan derivasyonlarda), 0.1 mV birimiyle.
  - Önceden ST depresyonu olanlarda "daha da çöken" ST.

## 5. Sonlanımlar
- **Birincil:** 30 gün içinde tüm nedenlere bağlı ölüm (indeks zamanından itibaren).
- **İkincil:**
  - 1 yıllık mortalite (Cox modeli).
  - Miyokard hasarı: indeksten −24 ile +48 saat arasında troponin en yüksek değeri > URL ve en düşük değere göre ≥ %20 artış. Yalnız troponini ölçülenlerde değerlendirilir. Ölçülme olasılığına göre ters olasılık ağırlıklandırması (IPW) uygulanır; ölçülme modeli tüm kovaryatları içerir.

## 6. Kovaryatlar (önceden tanımlı düzeltme seti)
- **Demografi:** yaş, cinsiyet, indeks yılı.
- **Anemi:** son Hb (spline), ΔHb, baz Hb, kanama kodu (D62 veya majör kanama; v3 07 tanımları).
- **Hastalık ağırlığı (indeksten ±12 saat içindeki en yakın değer):** kreatinin, BUN, Na, K, albümin, lökosit, trombosit, bilirubin, INR, laktat, glukoz, bikarbonat.
  - Eksik veriler çoklu atamayla (m=10) doldurulur.
  - Laktat ve INR için ayrıca "ölçülmüş mü" göstergesi eklenir; ölçülmüş olmaları bilgi taşıyabilir.
- **Kardiyak:**
  - KAH, kalp yetersizliği, KBH, diyabet ve Elixhauser skoru. MGH'de indeksten önce tarihlenmiş ICD kodları kullanılır. MIMIC'te yalnız önceki yatışların kodları kullanılır, çünkü aynı yatışın kodları sonlanımla iç içe geçer.
  - NT-proBNP ölçülmüş olma göstergesi.
- **EKG:** baz ST60_min, baz QRS süresi, son kalp hızı (bkz. §9 karar 1).

## 7. Analiz
1. **Birincil:** Lojistik regresyon (30 günlük ölüm) ~ yeni ST depresyonu + §6 kovaryatları. Sonuç: düzeltilmiş OR ve 95% CI.
   - **MGH başarı ölçütü:** Düzeltilmiş OR'nin 95% CI alt sınırı > 1.
   - **MIMIC ölçütü:** §9 karar 2.
2. **Özgüllük (anemi intoleransı iddiası için zorunlu):** Aynı model kontrast kohortunda kurulur. Veriler birleştirilip "yeni ST depresyonu × anemik kohort" etkileşimi test edilir; sonuç OR oranı ve 95% CI olarak verilir.
   - Etkileşim > 1 ise "anemiye özgü intolerans" yorumu desteklenir.
   - Etkileşim ≈ 1 ise sonuç "ST depresyonu genel bir risk göstergesidir; anemiye özgü değildir" olarak raporlanır.
3. **Ek değer:** Ağırlık modeline (§6, ST'siz) yeni ST depresyonu eklendiğinde ΔC-istatistiği (bootstrap 2000), kalibrasyon ve sürekli NRI (ikincil).
4. **Önceden tanımlı alt gruplar:**
   - son Hb ≤7 ile 7–8;
   - KAH var/yok;
   - kanama kodu var/yok;
   - kalp hızı artışı ≥10 vuru/dk var/yok.
5. **Duyarlılık analizleri:**
   - kalp hızı düzeltmesi yok/var (§9 karar 1'in tersi);
   - tam veri (complete-case);
   - yalnız 500 Hz kayıtlar;
   - 1 mm eşiği;
   - ilk çift yerine rastgele seçilen bir uygun çift.
6. **Derin öğrenme (ikincil, ana iddiayı taşımaz):**
   - Çift-EKG modeli MGH'nin %70 hasta-ayrık eğitim bölümünde eğitilir. Eğitimde v3 çift kohortunun tamamı kullanılır; indeks hastalar eğitimden çıkarılır.
   - Hedef: 30 günlük ölüm.
   - Değerlendirme MGH'nin kalan %30'unda ve MIMIC'te yapılır. Sorulan soru: §6 modeline ΔC katkısı.

## 8. Süreç kuralları
- **Kod dondurma ve MIMIC:**
  - Kod, MGH sonuç analizinden önce dondurulur ve hash'leri §11'e yazılır.
  - MIMIC inşa betiği yalnız MIMIC'e zorunlu farkları içerir ve bu farklar günlüğe yazılır.
  - MIMIC analizi, donmuş kodla tek seferde çalıştırılır.
- **Veri kullanım sözleşmeleri:** BDSP ve PhysioNet şartlarının teyidi (v3 §11.8) bu çalışmayı da kapsar.
- **Kısıtlar (makalede yazılacak):**
  - MGH'de yaşam bulguları (TA, SpO2) ve yoğun bakım bilgisi yok.
  - EKG ve troponin endikasyonla istenmiş.
  - Transfüzyon verisi yok; bu nedenle tedavi etkisi iddiası yapılmaz.
  - MIMIC'te ölüm tarihi sansürlü.

## 9. Kilitlenmiş kararlar (2026-10-02, kullanıcı onayladı)
1. **Son kalp hızı birincil modelde olsun mu?** Kalp hızı hem ağırlık göstergesi hem de intoleransın aracısı. Öneri: birincil modele **dahil edilsin** (muhafazakâr yaklaşım: "ST, taşikardinin ötesinde bilgi veriyor mu?"). Duyarlılık analizinde çıkarılır.
2. **MIMIC replikasyon ölçütü:** MIMIC'te ~350 ölüm ve ~63 maruz ölüm var. %80 güçle saptanabilecek en küçük OR yaklaşık 1.6. Öneri: replikasyon = nokta tahmini > 1 **ve** MGH ile heterojenlik p > 0.05. Birleşik tahmin birincil raporlanır. MIMIC'in tek başına anlamlılığı şart koşulmaz ve bu güç kısıtı makalede açıkça yazılır.
3. **Özgüllük testi zorunlu mu?** Öneri: **evet, ana iddianın parçası** olsun. "Anemi intoleransı" demek için etkileşim > 1 olmalı. Aksi halde makale yalnızca "akut anemide ST depresyonu kötü prognoz göstergesidir" der.

## 10. Değişiklik günlüğü
- 2026-10-02: Taslak v0.1 oluşturuldu. Kovaryat laboratuvarlarının (laktat, lökosit, trombosit, bilirubin, INR, NT-proBNP, glukoz, bikarbonat) iki merkezde de mevcut olduğu doğrulandı.
- 2026-10-02: **v1.0 KİLİTLENDİ.** Kararlar: (1) son kalp hızı birincil modelde yer alır; (2) MIMIC replikasyonu = nokta tahmini > 1 ve heterojenlik p > 0.05, birincil rapor birleşik tahmin; (3) özgüllük (etkileşim) testi ana iddianın parçası.
- 2026-10-02 **Kod donduruldu** (MGH sonuç analizinden önce): `heedb/10_build_tolerance_mgh.py` sha1 df7fae37, `heedb/11_tolerance_analysis.py` sha1 5954a4ad.
  - Kod, sonuç değişkenleri karıştırılmış sahte veride duman testinden geçirildi. Gerçek sonuçlarla model çalıştırılmadı.
  - **Sonuç görülmeden yapılan zorunlu değişiklikler:**
    - (a) Baz Hb kovaryat setinden çıkarıldı. ΔHb = baz Hb − son Hb tanımı gereği, son Hb spline'ı ile birlikte tam doğrusal bağımlılık oluşturuyordu; bilgi kaybı yok.
    - (b) Sürekli kovaryatlar sayısal kararlılık için standartlaştırıldı; maruziyetin OR'si etkilenmez.
    - (c) Analiz kümesinde 10'dan az pozitif veya negatif değeri olan ikili kovaryatlar modelden çıkarılır (sonuçtan bağımsız kural).
    - (d) Ölüm tarihi indeks EKG'den önce görünen 9 hasta çıkarıldı.
    - (e) Troponin "ölçülmüş" tanımı, artış değerlendirilebilsin diye ≥2 değer olarak uygulandı.
    - (f) NT-proBNP "ölçülmüş" göstergesi için ±48 saat penceresi kullanıldı.
  - **Ayrıca bildirilmesi gereken:** Kohort inşa kontrolünde karşılaştırma kohortunun ham 30 günlük ölüm oranı (%6.7) görüldü. Herhangi bir model çalıştırılmadı.
  - **Derin öğrenme (ikincil):** Ayrı bir adım olarak çalıştırılacak.
- 2026-10-02: MIMIC inşa betiği `heedb/10b_build_tolerance_mimic.py` (sha1 3bfac8b4) yazıldı ve çalıştırıldı: 1.928 anemik indeks hasta, 7.810 karşılaştırma hastası. Sonuç analizi henüz yapılmadı.
  - MIMIC'e özgü farklar:
    - komorbiditeler yalnız önceki yatışların kodlarından;
    - kanama kodu aynı yatış düzeyinde (D62 %33.5; MGH'de ±3 gün tanımıyla %0.2);
    - takvim yılı anchor_year_group orta noktası üzerinden yaklaşık;
    - troponin T (51003);
    - 30 günlük takip tam kabul edildi.
- 2026-10-02 **SONUÇLAR** (donmuş kod 5954a4ad, iki merkezde aynı):
  - **Birincil (MGH):** OR 1.35 [1.16–1.57], 17.977 hasta / 3.052 ölüm. **Başarı ölçütü sağlandı.**
  - **MIMIC:** OR 1.27 [0.86–1.87], heterojenlik p=0.77. **Önceden tanımlı replikasyon ölçütü sağlandı** (nokta tahmini > 1 ve heterojenlik yok).
  - **Birleşik tahmin:** OR 1.34 [1.16–1.54].
  - **Özgüllük:** Anemik/karşılaştırma OR oranı MGH'de 1.16 [0.95–1.42], MIMIC'te 1.32 [0.81–2.16]; birleşik 1.18 [0.98–1.42], p=0.076. **Anemiye özgü intolerans iddiası DESTEKLENMEDİ** (§9.3 kuralı). Sonuç "akut anemide yeni ST depresyonu kötü prognoz göstergesidir" olarak raporlanır.
  - **Ek değer:** ΔC MGH'de 0.0008 [0.0001–0.0017], MIMIC'te 0.0004. NRI anlamsız. Klinik olarak ihmal edilebilir.
  - **Doz-yanıt yok:** 1 mm eşiğiyle OR 1.02 (MGH) ve 0.74 (MIMIC). Sürekli ST yükü OR 1.02.
  - **Kalp hızı aracılığı:** HR düzeltmesi yapılmayınca OR 1.51 (MGH) ve 1.58 (MIMIC).
  - **Alt gruplar:** Kanama alt grubunda etkileşim iki merkezde de tersine (OR oranı 0.28 ve 0.35, ikisi de p<0.05); kanayanlarda ST depresyonu mortaliteyle ilişkili değil. Bu 8 alt grup testinden biri; keşifsel olarak yorumlanmalı.
  - **1 yıllık mortalite:** HR 1.16 [1.06–1.27] (MGH), 1.10 [0.88–1.38] (MIMIC).
  - **Miyokard hasarı (IPW):** MGH OR 1.28 [0.96–1.71]; MIMIC'te tahmin edilemedi (küçük n).
