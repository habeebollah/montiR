# montiR
montiR (Model biOmass diNamik singkaT dI R) merupakan kumpulan model biomas dinamik dengan pendekatan time series fitting. Tool pemodelan yang tersedia disini antara lain:

(1) Surplus produksi dengan asumsi equilibrium

Tool ini menghitung jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) dengan asumsi equilibrium untuk model Schaefer dan Fox. Pendekatan ini ditampilkan disini hanya untuk tujuan edukasi sebagai contoh model yang akan memberikan estimasi MSY dan Emsy yang lebih tinggi, sehingga sangat tidak disarankan untuk dijadikan sebagai panduan dalam pengambilan kebijakan perikanan. Overestimasi reference point pada kondisi ketika status perikanan sedang dalam kondisi overexploited akan memberikan ilusi bahwa stok ikan masih banyak, sehingga dapat merugikan pelaku perikanan karena jumlah tangkapan yang rendah dan merugikan stok ikan karena semakin tingginya pemanfaatan. Banyak dari kita yang masih menggunakan metode ini, meskipun sudah tidak disarankan untuk digunakan sejak 1980an.

(2) Surplus produksi dengan asumsi non-equilibrium menggunakan multiple regression

Metode multiple regression masih banyak digunakan untuk menghitung jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) dengan asumsi non-equilibrium untuk model Schaefer. Metode yang paling banyak digunakan disebut menghasilkan bias terhadap parameter yang diestimasi, sehingga menghasilkan perhitungan MSY dan Emsy yang juga bias. Revisi dari metode ini kurang banyak digunakan, selain itu semua metode yang menggunakan multiple regression disebut memiliki kelemahan jika menggunakan data yang memiliki tipe one way trip.  

(3) Surplus produksi dengan asumsi non-equilibrium menggunakan data fitting

(3.a) Data plotting

Langkah paling penting sebelum melakukan analisis data adalah memeriksa apakah data yang akan digunakan memenuhi persyaratan dan asumsi yang dibutuhkan untuk analisis biomass dynamic model, termasuk memilih jenis langkah apa yang harus dilakukan ketika data yang dibutuhkan tidak memenuhi asumsi

(3.b) Surplus produksi dengan asumsi non-equilibrium menggunakan time series fitting

Tool ini melakukan estimasi parameter K, B0, r, q dan menentukan jumlah tangkapan ikan lestari (MSY), biomassa ikan lestari (Bmsy), serta upaya penangkapan ikan lestari (Emsy) menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox. Tool ini sudah disesuaikan untuk kebutuhan data yang terbatas (dapat mengakomodasi hilangnya input data upaya penangkapan) serta sudah memperhitungkan kesalahan dalam pengambilan data (observation error). Metode time series fitting disebut sebagai metode yang lebih baik dibandingkan dengan dua metode lain (metode equilibrium dan multiple regression) yang digunakan untuk melakukan estimasi parameter dalam model surplus produksi. Sebagian kecil dari kita sudah menggunakan metode ini, tetapi pelaksanaan analisisnya masih perlu disempurnakan agar tidak berakibat pada kurang tepatnya perhitungan MSY, Bmsy dan Emsy

(3.c) Menghitung reference point untuk pengelolaan

Tool ini menghasilkan jumlah stok ikan yang lestari (Bmsy), jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) untuk model Schaefer dan Fox.

(3.d) Menghitung standard error dari reference point

Tool ini mengestimasi jumlah stok ikan yang lestari (Bmsy), jumlah tangkapan ikan lestari (MSY) dan upaya penangkapan ikan lestari (Emsy) serta menghitung standard error menggunakan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox.

(4) Mengambil data prior untuk parameter pertumbuhan r sebagai input pendugaan stok di tingkat spesies menggunakan bayesian

Data prior setiap spesies untuk parameter pertumbuhan r diambil dari database fishbase dan sealifebase untuk mendukung pendugaan stok di tingkat spesies. Contoh langkah lanjutan untuk melakukan estimasi parameter surplus production menggunakan pendekatan bayesian dan data runut waktu dengan asumsi non-equilibrium untuk model Schaefer dan Fox sehingga menghasilkan pendugaan stok yang lebih akurat di tingkat spesies juga diberikan.

(5) Membuat proyeksi atas kebijakan reference point berdasar tingkat pemanfaatan perikanan

Tool ini akan membuat grafik proyeksi biomass per biomass at msy (B/Bmsy) dan fishing per fishing at msy (F/Fmsy) sebagai panduan untuk melihat kebijakan yang akan dibuat berdasarkan Bmsy, MSY dan Emsy sebagai reference point. Proyeksi dibuat dengan pendekatan deterministic secara default, dan terdapat opsi untuk tujuan stochastic.

## Installation

Anda dapat mengunduh dan menginstal package ini melalui link berikut:

``` r
devtools::install_github("habeebollah/montiR")
```

## Contoh

Berikut adalah contoh sederhana penggunaan montiR

```{r example}
library(montiR)
## membuat input data
df <- data.frame(tahun=c(...),
                 tangkapan=c(...),
                 upaya=c(...))
                 
## melihat bentuk data
plotInit(df=df)

## mencari pilihan parameter awal sebagai input untuk estimasi
K <- 1000
B0 <- K
r <- 0.3
q <- 0.00025

inpars <- c(log(K), log(B0), log(r), log(q))
SPparS(inpars=inpars, df=goodcontrast)
```

