# FALLTRACE

**Düşen uzay enkazını dünyanın deprem sensörleriyle izle.**

Bir uydu ya da roket kademesi atmosfere geri girdiğinde, oluşturduğu ses patlaması yeri
sarsar. Binlerce açık sismometre bunu kaydeder, ama bu veriye kimse bakmaz. FALLTRACE bu
kayıtlardan cismin yörüngesini, yönünü, hızını ve enkazın düşebileceği koridoru çıkarır.
Her sonucu dürüst bir hata payıyla verir.

> Durum: **erken geliştirme (pre-alpha)**. Fizik çekirdeği ve doğrulama testleri hazır.
> Gerçek olaylarla doğrulama bir sonraki adım.

## Neden önemli?

- Her gün birkaç Starlink uydusu atmosfere giriyor, roket kademeleri de düzenli olarak düşüyor.
- Kontrolsüz bir düşüşün zamanı çoğu zaman ancak saatler mertebesinde tahmin edilebiliyor.
  Cisim yaklaşık 7,5 km/s hızla gittiği için bir saatlik belirsizlik, yer izinde on binlerce
  kilometre demek.
- Kanıt zaten kaydediliyor. 2025'te bir Falcon 9 üst kademesi Polonya'ya düştü. 2024'te
  Shenzhou-15 modülünün Kaliforniya üzerinden geçişi 100'den fazla sismometreye kaydedildi.

## Hızlı başlangıç

```bash
pip install -e ".[all,dev]"
falltrace demo
falltrace fit examples/synthetic_picks.csv -o iz.geojson
```

## Dürüstlük ilkesi

1. Her sayının bir belirsizliği var. Testler, bu belirsizliklerin gerçekçi olup olmadığını kontrol ediyor.
2. Bilinen sapmalar ölçülüp yayınlanıyor. Örneğin ses hızını %5 yanlış varsaymak, konumu iz
   boyunca yaklaşık 100 km kaydırıyor. Bu yüzden ses hızı varsayılan olarak veriden tahmin ediliyor.
3. Başarısız sonuçlar da [doğrulama sayfasında](docs/validation.md) paylaşılıyor.
4. Sentetik veride başarılı olmak, gerçek dünyada doğrulanmış olmak demek değil.

Detaylar için İngilizce [README](README.md), [yöntem](docs/method.md) ve
[doğrulama](docs/validation.md) sayfalarına bakabilirsin.
