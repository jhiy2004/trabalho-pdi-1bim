#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QMessageBox>
#include <QMouseEvent>

struct DCTTask {
    const QImage& img;
    int m, n;
    std::vector<std::vector<float>>& dct;
    const std::vector<std::vector<float>>& cos_k_i;
    const std::vector<std::vector<float>>& cos_l_j;

    DCTTask(const QImage& img_,
            int m_, int n_,
            std::vector<std::vector<float>>& dct_,
            const std::vector<std::vector<float>>& cos_k_i_,
            const std::vector<std::vector<float>>& cos_l_j_)
        : img(img_), m(m_), n(n_), dct(dct_),
        cos_k_i(cos_k_i_), cos_l_j(cos_l_j_) {}

    void operator()(int i) const {
        float ci = (i == 0) ? 1.0f / sqrt(m) : sqrt(2.0f / m);

        for (int j = 0; j < n; ++j) {
            float cj = (j == 0) ? 1.0f / sqrt(n) : sqrt(2.0f / n);
            float sum = 0.0f;

            for (int k = 0; k < m; ++k) {
                float cos_ki = cos_k_i[k][i];

                for (int l = 0; l < n; ++l) {
                    float gray = qGray(img.pixel(l, k));
                    float cos_lj = cos_l_j[l][j];
                    sum += gray * cos_ki * cos_lj;
                }
            }

            dct[i][j] = ci * cj * sum;
        }
    }
};

struct IDCTTask {
    const std::vector<std::vector<float>>& dct;
    int m, n;
    QImage* result;
    const std::vector<std::vector<float>>& cos_i_k;
    const std::vector<std::vector<float>>& cos_j_l;

    IDCTTask(
        const std::vector<std::vector<float>>& dct_,
        QImage* result_,
        int m_, int n_,
        const std::vector<std::vector<float>>& cos_i_k_,
        const std::vector<std::vector<float>>& cos_j_l_)
        : dct(dct_), result(result_), m(m_), n(n_),
        cos_i_k(cos_i_k_), cos_j_l(cos_j_l_) {}

    void operator()(int i) const {
        for (int j = 0; j < n; ++j) {
            float sum = 0.0f;

            for (int k = 0; k < m; ++k) {
                float ck = (k == 0) ? 1.0f / sqrt(m) : sqrt(2.0f / m);
                float cos_ik = cos_i_k[i][k];

                for (int l = 0; l < n; ++l) {
                    float cl = (l == 0) ? 1.0f / sqrt(n) : sqrt(2.0f / n);
                    float cos_jl = cos_j_l[j][l];
                    sum += ck * cl * dct[k][l] * cos_ik * cos_jl;
                }
            }

            int gray = qBound(0, static_cast<int>(sum + 0.5f), 255);
            QRgb* line = reinterpret_cast<QRgb*>(result->scanLine(i));
            line[j] = qRgb(gray, gray, gray);
        }
    }
};


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->output_image->setMouseTracking(true);

    ui->output_image->installEventFilter(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

QImage img;
QImage res;

std::vector<std::vector<float>> dct_output;
std::vector<std::vector<float>> dct_input;

void MainWindow::on_actionEscala_Cinza_triggered()
{
    sobelAtivo = false;
    QImage mod = img;

    int w = ui->output_image->width();
    int h = ui->output_image->height();

    for(int i=0; i < mod.height(); i++){
        for(int j=0; j < mod.width(); j++){
            QRgb color = mod.pixel(j, i);
            int res = static_cast<int>(0.299f * qRed(color) + 0.587f * qGreen(color) + 0.114f * qBlue(color));

            QRgb cinza = qRgb(res, res, res);
            mod.setPixel(j, i, cinza);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);

    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}

void MainWindow::on_actionCarregar_imagem_triggered()
{
    sobelAtivo = false;
    QString filename = QFileDialog::getOpenFileName(this, "Open a file", ".");
    img.load(filename);

    int w = ui->input_image->width();
    int h = ui->input_image->height();
    QPixmap pix(filename);
    ui->input_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
}

void MainWindow::on_actionEscala_Cinza_Inversa_triggered()
{
    sobelAtivo = false;
    QImage mod = img;

    int w = ui->output_image->width();
    int h = ui->output_image->height();

    for(int i=0; i < mod.height(); i++){
        for(int j=0; j < mod.width(); j++){
            QRgb color = mod.pixel(j, i);
            int res = 255 - static_cast<int>(0.299f * qRed(color) + 0.587f * qGreen(color) + 0.114f * qBlue(color));

            QRgb cinza = qRgb(res, res, res);
            mod.setPixel(j, i, cinza);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);

    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionColorida_Inversa_triggered()
{
    sobelAtivo = false;
    QImage mod = img;

    int w = ui->output_image->width();
    int h = ui->output_image->height();

    for(int i=0; i < mod.height(); i++){
        for(int j=0; j < mod.width(); j++){
            QRgb color = mod.pixel(j, i);
            QRgb inv_color = qRgb(255 - qRed(color), 255 - qGreen(color), 255 - qBlue(color));
            mod.setPixel(j, i, inv_color);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);

    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionSeparar_canais_de_cores_triggered()
{
    sobelAtivo = false;
    QImage red_channel = img;
    QImage green_channel = img;
    QImage blue_channel = img;


    for(int i=0; i < img.height(); i++){
        for(int j=0; j < img.width(); j++){
            QRgb color = img.pixel(j, i);

            red_channel.setPixel(j, i, qRgb(qRed(color), 0, 0));
            green_channel.setPixel(j, i, qRgb(0, qGreen(color), 0));
            blue_channel.setPixel(j, i, qRgb(0, 0, qBlue(color)));
        }
    }

    QPixmap red_pix = QPixmap::fromImage(red_channel);
    QPixmap green_pix = QPixmap::fromImage(green_channel);
    QPixmap blue_pix = QPixmap::fromImage(blue_channel);

    int w = ui->channel_1->width();
    int h = ui->channel_1->height();
    ui->channel_1->setPixmap(red_pix.scaled(w, h, Qt::KeepAspectRatio));

    w = ui->channel_2->width();
    h = ui->channel_2->height();
    ui->channel_2->setPixmap(green_pix.scaled(w, h, Qt::KeepAspectRatio));

    w = ui->channel_3->width();
    h = ui->channel_3->height();
    ui->channel_3->setPixmap(blue_pix.scaled(w, h, Qt::KeepAspectRatio));
}


void MainWindow::on_actionSal_e_pimenta_triggered()
{
    sobelAtivo = false;
    if (img.isNull()) {
        qDebug() << "Erro: Imagem não carregada corretamente!";
        return;
    }

    QImage mod = img;

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    int w = ui->output_image->width();
    int h = ui->output_image->height();

    int tamanho = img.width() * img.height();
    int nRuido = tamanho / 10; // 10% dos pixels serão ruídos
    int x, y, c;

    for(int i = 0; i < nRuido; i++) {
        x = QRandomGenerator::global()->bounded(img.width());
        y = QRandomGenerator::global()->bounded(img.height());
        c = QRandomGenerator::global()->bounded(2) * 255; // 0 (pimenta) ou 255 (sal)

        QRgb new_c = qRgb(c, c, c);
        mod.setPixel(x, y, new_c);
    }

    QPixmap pix = QPixmap::fromImage(mod);
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}

void MainWindow::on_output_to_input_btn_clicked()
{
    sobelAtivo = false;
    QPixmap output_pix = ui->output_image->pixmap();
    if (!output_pix || output_pix.isNull())
        return;

    img = res;
    if(!dct_output.empty()){
        dct_input = dct_output;
    }
    QPixmap pix = QPixmap::fromImage(img);

    int w = ui->input_image->width();
    int h = ui->input_image->height();
    ui->input_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
}



void MainWindow::on_actionFiltro_mediana_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int matrix[9];
            matrix[0] = qRed(img.pixel(j-1, i-1));
            matrix[1] = qRed(img.pixel(j-1, i));
            matrix[2] = qRed(img.pixel(j-1, i+1));
            matrix[3] = qRed(img.pixel(j, i-1));
            matrix[4] = qRed(img.pixel(j, i));
            matrix[5] = qRed(img.pixel(j, i+1));
            matrix[6] = qRed(img.pixel(j+1, i-1));
            matrix[7] = qRed(img.pixel(j+1, i));
            matrix[8] = qRed(img.pixel(j+1, i+1));

            std::sort(matrix, matrix + 9);

            int median = matrix[4];

            QRgb color = qRgb(median, median, median);
            mod.setPixel(j-1, i-1, color);
        }
    }
    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}

void MainWindow::on_actionFiltro_m_dia_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int sum=0;
            sum += qRed(img.pixel(j-1, i-1));
            sum += qRed(img.pixel(j-1, i));
            sum += qRed(img.pixel(j-1, i+1));
            sum += qRed(img.pixel(j, i-1));
            sum += qRed(img.pixel(j, i));
            sum += qRed(img.pixel(j, i+1));
            sum += qRed(img.pixel(j+1, i-1));
            sum += qRed(img.pixel(j+1, i));
            sum += qRed(img.pixel(j+1, i+1));

            int media = sum/9;

            QRgb color = qRgb(media, media, media);
            mod.setPixel(j-1, i-1, color);
        }
    }
    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionBinarizar_triggered()
{
    sobelAtivo = false;
    bool ok=false;
    int threshold = QInputDialog::getInt(this, tr("Digite o limiar"),
                                tr("Limiar:"), 128, 0, 255, 1, &ok);

    if(ok){
        QImage mod = img;

        for(int i=0; i < img.height(); ++i){
            for(int j=0; j < img.width(); ++j){
                int color = qRed(img.pixel(j, i));
                color = (color >= threshold) ? 255 : 0;

                mod.setPixel(j, i, qRgb(color, color, color));
            }
        }
        QPixmap pix = QPixmap::fromImage(mod);

        int w = ui->output_image->width();
        int h = ui->output_image->height();
        ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
        res = mod;
    }
}

void MainWindow::on_actionLaplaciano_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int sum = 0;
            sum += -qGray(img.pixel(j-1, i));
            sum += -qGray(img.pixel(j, i-1));
            sum +=  4*qGray(img.pixel(j, i));
            sum += -qGray(img.pixel(j, i+1));
            sum += -qGray(img.pixel(j+1, i));

            int new_color = qBound(0, sum, 255);

            mod.setPixel(j-1, i-1, qRgb(new_color, new_color, new_color));
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionEqualiza_o_triggered()
{
    sobelAtivo = false;
    int largura = img.width();
    int altura = img.height();
    int tamanho = largura * altura;

    QVector<int> histograma(256, 0);
    for (int y = 0; y < altura; ++y) {
        for (int x = 0; x < largura; ++x) {
            int intensidade = qGray(img.pixel(x, y));
            histograma[intensidade]++;
        }
    }

    QVector<int> cdf(256, 0);
    cdf[0] = histograma[0];
    for (int i = 1; i < 256; ++i) {
        cdf[i] = cdf[i - 1] + histograma[i];
    }

    QVector<int> novaIntensidade(256, 0);
    for (int i = 0; i < 256; ++i) {
        novaIntensidade[i] = qRound(((cdf[i] - cdf[0]) / double(tamanho - cdf[0])) * 255);
    }

    QImage imagemEqualizada(largura, altura, QImage::Format_Grayscale8);
    for (int y = 0; y < altura; ++y) {
        for (int x = 0; x < largura; ++x) {
            int intensidade = qGray(img.pixel(x, y));
            int novaInt = novaIntensidade[intensidade];
            imagemEqualizada.setPixel(x, y, qRgb(novaInt, novaInt, novaInt));
        }
    }

    QPixmap pix = QPixmap::fromImage(imagemEqualizada);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = imagemEqualizada;
}


void MainWindow::on_actionCanal_Vermelho_em_Cinza_triggered()
{
    sobelAtivo = false;
    int largura = img.width();
    int altura = img.height();
    QImage red_channel = img;
    QImage green_channel = img;
    QImage blue_channel = img;

    for(int i = 0; i < altura; ++i){
        for(int j = 0; j < largura; ++j){
            QRgb color = img.pixel(j, i);

            red_channel.setPixel(j, i, qRgb(qRed(color), qRed(color), qRed(color)));
            green_channel.setPixel(j, i, qRgb(qGreen(color), qGreen(color), qGreen(color)));
            blue_channel.setPixel(j, i, qRgb(qBlue(color), qBlue(color), qBlue(color)));
        }
    }

    QPixmap red_pix = QPixmap::fromImage(red_channel);
    QPixmap green_pix = QPixmap::fromImage(green_channel);
    QPixmap blue_pix = QPixmap::fromImage(blue_channel);

    int l = ui->channel_1->width();
    int a = ui->channel_1->height();
    ui->channel_1->setPixmap(red_pix.scaled(l, a, Qt::KeepAspectRatio));

    l = ui->channel_2->width();
    a = ui->channel_2->height();
    ui->channel_2->setPixmap(green_pix.scaled(l, a, Qt::KeepAspectRatio));

    l = ui->channel_3->width();
    a = ui->channel_3->height();
    ui->channel_3->setPixmap(blue_pix.scaled(l, a, Qt::KeepAspectRatio));
}


void MainWindow::on_actionRGB_para_HSV_triggered()
{
    sobelAtivo = false;
    bool ok;
    int R = QInputDialog::getInt(this, "Entrada RGB", "Digite o valor de R (0-255):", 0, 0, 255, 1, &ok);
    if (!ok) return;
    int G = QInputDialog::getInt(this, "Entrada RGB", "Digite o valor de G (0-255):", 0, 0, 255, 1, &ok);
    if (!ok) return;
    int B = QInputDialog::getInt(this, "Entrada RGB", "Digite o valor de B (0-255):", 0, 0, 255, 1, &ok);
    if (!ok) return;

    // Normalizando RGB
    float r = R / 255.0;
    float g = G / 255.0;
    float b = B / 255.0;

    // Obtendo V
    float V = std::max({r, g, b});

    // Obtendo S
    float min_val = std::min({r, g, b});
    float delta = V - min_val;

    float S;

    if(V == 0){
        S = 0;
    }else{
        S = delta / V;
    }

    // Obtendo H
    float H;

    if(delta == 0){
        H = 0;
    }else{
        if(V == r){
            H = 60 * fmod(((g - b) / delta), 6);
        }
        if(V == g){
            H = 60 * ((b - r) / delta + 2);
        }
        if(V == b){
            H = 60 * ((r - g) / delta + 4);
        }
    }

    if(H < 0){
        H += 360;    sobelAtivo = false;
        QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

        // Garante que o formato da imagem é compatível com setPixel()
        if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
            mod = mod.convertToFormat(QImage::Format_RGB32);
        }

        for(int i=1; i < img.height() - 1; ++i){
            for(int j=1; j < img.width() - 1; ++j){
                int sum = 0;
                sum += -qGray(img.pixel(j-1, i));
                sum += -qGray(img.pixel(j, i-1));
                sum +=  4*qGray(img.pixel(j, i));
                sum += -qGray(img.pixel(j, i+1));
                sum += -qGray(img.pixel(j+1, i));

                int new_color = qBound(0, sum, 255);

                mod.setPixel(j-1, i-1, qRgb(new_color, new_color, new_color));
            }
        }

        QPixmap pix = QPixmap::fromImage(mod);

        int w = ui->output_image->width();
        int h = ui->output_image->height();
        ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
        res = mod;
    }

    // Transformo em porcentagem
    S = S * 100;
    V = V * 100;

    QString result = QString("HSV:\nH: %1\nS: %2\nV: %3").arg(H).arg(S).arg(V);
    QMessageBox::information(this, "Resultado HSV", result);
}


void MainWindow::on_actionHSV_para_RGB_triggered()
{
    sobelAtivo = false;
    bool ok;
    int H = QInputDialog::getInt(this, "Entrada HSV", "Digite o valor de H (0-360):", 0, 0, 360, 1, &ok);
    if (!ok) return;
    int Si = QInputDialog::getInt(this, "Entrada HSV", "Digite o valor de S (0-100):", 0, 0, 100, 1, &ok);
    if (!ok) return;
    int Vi = QInputDialog::getInt(this, "Entrada HSV", "Digite o valor de V (0-100):", 0, 0, 100, 1, &ok);
    if (!ok) return;

    float S = Si / 100.0;
    float V = Vi / 100.0;

    float C = V * S;
    float X = C * (1 - fabs(fmod((H / 60.0), 2) - 1));
    float m = V - C;

    float r = 0, g = 0, b = 0;

    if (0 <= H && H < 60) {
        r = C;
        g = X;
        b = 0;
    } else if (60 <= H && H < 120) {
        r = X;
        g = C;
        b = 0;
    } else if (120 <= H && H < 180) {
        r = 0;
        g = C;
        b = X;
    } else if (180 <= H && H < 240) {
        r = 0;
        g = X;
        b = C;
    } else if (240 <= H && H < 300) {
        r = X;
        g = 0;
        b = C;
    } else if (300 <= H && H < 360) {
        r = C;
        g = 0;
        b = X;
    }

    int R = qRound((r + m) * 255);
    int G = qRound((g + m) * 255);
    int B = qRound((b + m) * 255);

    QString result = QString("RGB:\nR: %1\nG: %2\nB: %3").arg(R).arg(G).arg(B);
    QMessageBox::information(this, "Resultado RGB", result);
}


void MainWindow::on_actionCompress_o_de_Escala_Din_mica_triggered()
{
    sobelAtivo = false;
    bool ok;
    double c = QInputDialog::getDouble(this, "Entrada", "Digite o valor de C (0-1):", 1.0, 0.0, 1.0, 4, &ok);
    if (!ok) return;
    double y = QInputDialog::getDouble(this, "Entrada", "Digite o valor de Y:", 1.0, -1000000.0, 1000000.0, 4, &ok);
    if (!ok) return;

    int largura = img.width();
    int altura = img.height();
    QImage image = img;

    for(int i = 0; i < altura; ++i){
        for(int j = 0; j < largura; ++j){
            int intensidade = qGray(img.pixel(j, i));
            double r = intensidade / 255.0;
            double S = c * pow(r, y);
            int valorFinal = qBound(0, static_cast<int>(S * 255.0), 255);
            image.setPixel(j, i, qRgb(valorFinal, valorFinal, valorFinal));
        }
    }

    QPixmap pix = QPixmap::fromImage(image);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = image;
}


void MainWindow::on_actionBordas_por_Sobel_triggered()
{
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    magnitudes.clear();
    angles.clear();

    int min_magnitude = INT_MAX;
    int max_magnitude = INT_MIN;

    for (int i = 1; i < img.height() - 1; ++i) {
        std::vector<int> row;
        std::vector<float> angleRow;

        for (int j = 1; j < img.width() - 1; ++j) {
            int sumY = 0;
            sumY += -2 * qGray(img.pixel(j, i - 1));
            sumY += -1 * qGray(img.pixel(j - 1, i - 1));
            sumY += -1 * qGray(img.pixel(j + 1, i - 1));
            sumY += 2 * qGray(img.pixel(j, i + 1));
            sumY += qGray(img.pixel(j + 1, i + 1));
            sumY += qGray(img.pixel(j - 1, i + 1));

            int sumX = 0;
            sumX += -2 * qGray(img.pixel(j - 1, i));
            sumX += -1 * qGray(img.pixel(j - 1, i - 1));
            sumX += -1 * qGray(img.pixel(j - 1, i + 1));
            sumX += 2 * qGray(img.pixel(j + 1, i));
            sumX += qGray(img.pixel(j + 1, i + 1));
            sumX += qGray(img.pixel(j + 1, i - 1));

            int magnitude = qSqrt(sumX * sumX + sumY * sumY);

            float angle = qAtan2(sumY, sumX) * (180.0 / M_PI); // -180° a 180°
            if (angle < 0)
                angle += 360.0; // converte para 0° a 360°

            row.push_back(magnitude);
            angleRow.push_back(angle);

            if (magnitude < min_magnitude) min_magnitude = magnitude;
            if (magnitude > max_magnitude) max_magnitude = magnitude;
        }

        magnitudes.push_back(row);
        angles.push_back(angleRow);
    }

    for (int i = 0; i < mod.height(); ++i) {
        for (int j = 0; j < mod.width(); ++j) {
            int magnitude = magnitudes[i][j];
            int norm = 0;
            if (max_magnitude != min_magnitude) {
                norm = (magnitude - min_magnitude) * 255 / (max_magnitude - min_magnitude);
            }
            QRgb color = qRgb(norm, norm, norm);
            magnitudes[i][j] = norm;
            mod.setPixel(j, i, color);
        }
    }

    sobelImage = mod;
    sobelAtivo = true;

    QPixmap pix = QPixmap::fromImage(mod);
    ui->output_image->setPixmap(pix.scaled(ui->output_image->width(), ui->output_image->height(), Qt::KeepAspectRatio));
    res = mod;
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (obj == ui->output_image && sobelAtivo) {

        if (event->type() == QEvent::MouseMove) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);

            if (sobelImage.isNull()) return true;

            QPoint pos = mouseEvent->pos();

            QSize labelSize = ui->output_image->size();
            QSize imageSize = sobelImage.size();

            double scale = qMin(
                static_cast<double>(labelSize.width()) / imageSize.width(),
                static_cast<double>(labelSize.height()) / imageSize.height()
                );

            QSize scaledSize = imageSize * scale;

            int offsetX = (labelSize.width() - scaledSize.width()) / 2;
            int offsetY = (labelSize.height() - scaledSize.height()) / 2;

            if (pos.x() < offsetX || pos.x() >= offsetX + scaledSize.width() ||
                pos.y() < offsetY || pos.y() >= offsetY + scaledSize.height()) {
                ui->magnitude_label->clear();
                return true;
            }

            int imgX = (pos.x() - offsetX) * imageSize.width() / scaledSize.width();
            int imgY = (pos.y() - offsetY) * imageSize.height() / scaledSize.height();

            if (imgX >= 0 && imgX < imageSize.width() &&
                imgY >= 0 && imgY < imageSize.height()) {

                int mag = magnitudes[imgY][imgX];
                float angle = angles[imgY][imgX];
                ui->magnitude_label->setText(QString("Magnitude: %1 | Direção: %2°").arg(mag).arg(angle, 0, 'f', 1));
            } else {
                ui->magnitude_label->clear();
            }

            return true;
        }

        if (event->type() == QEvent::Leave) {
            ui->magnitude_label->clear();
            return true;
        }
    }

    return QMainWindow::eventFilter(obj, event);
}

void MainWindow::on_actionLimiariza_o_triggered()
{
    sobelAtivo = false;
    bool ok=false;
    int threshold = QInputDialog::getInt(this, tr("Digite o limiar"),
                                         tr("Limiar:"), 128, 0, 255, 1, &ok);

    if(ok){
        QImage mod = img;

        for(int i=0; i < img.height(); ++i){
            for(int j=0; j < img.width(); ++j){
                int color = qRed(img.pixel(j, i));
                color = (color >= threshold) ? color : 0;

                mod.setPixel(j, i, qRgb(color, color, color));
            }
        }
        QPixmap pix = QPixmap::fromImage(mod);

        int w = ui->output_image->width();
        int h = ui->output_image->height();
        ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
        res = mod;
    }
}

void MainWindow::on_actionDCT_triggered()
{
    const int m = img.height();
    const int n = img.width();

    dct_output.clear();
    dct_output.resize(m, std::vector<float>(n));

    // Precompute cosine tables
    std::vector<std::vector<float>> cos_k_i(m, std::vector<float>(m));
    std::vector<std::vector<float>> cos_l_j(n, std::vector<float>(n));

    for (int k = 0; k < m; ++k)
        for (int i = 0; i < m; ++i)
            cos_k_i[k][i] = std::cos((2 * k + 1) * i * M_PI / (2.0f * m));

    for (int l = 0; l < n; ++l)
        for (int j = 0; j < n; ++j)
            cos_l_j[l][j] = std::cos((2 * l + 1) * j * M_PI / (2.0f * n));

    DCTTask task(img, m, n, dct_output, cos_k_i, cos_l_j);

    QVector<int> rows(m);
    std::iota(rows.begin(), rows.end(), 0);

    QtConcurrent::blockingMap(rows, task);

    // Convert DCT result to grayscale image
    QImage mod(img.size(), QImage::Format_RGB32);
    for (int i = 0; i < m; ++i) {
        QRgb* line = reinterpret_cast<QRgb*>(mod.scanLine(i));
        for (int j = 0; j < n; ++j) {
            int val = std::clamp(static_cast<int>(dct_output[i][j]), 0, 255);
            line[j] = qRgb(val, val, val);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);
    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));

    res = mod;
}

void MainWindow::on_actionIDCT_triggered()
{
    const int m = img.height();
    const int n = img.width();

    QImage mod(img.size(), QImage::Format_RGB32);

    // Precompute cosine values
    std::vector<std::vector<float>> cos_i_k(m, std::vector<float>(m));
    std::vector<std::vector<float>> cos_j_l(n, std::vector<float>(n));

    for (int i = 0; i < m; ++i)
        for (int k = 0; k < m; ++k)
            cos_i_k[i][k] = std::cos((2 * i + 1) * k * M_PI / (2.0f * m));

    for (int j = 0; j < n; ++j)
        for (int l = 0; l < n; ++l)
            cos_j_l[j][l] = std::cos((2 * j + 1) * l * M_PI / (2.0f * n));

    // Prepare task and parallel map
    IDCTTask task(dct_input, &mod, m, n, cos_i_k, cos_j_l);
    QVector<int> rows(m);
    std::iota(rows.begin(), rows.end(), 0);
    QtConcurrent::blockingMap(rows, task);

    // Display result
    QPixmap pix = QPixmap::fromImage(mod);
    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionFiltragem_passa_baixa_DCT_triggered()
{
    bool ok=false;
    int cut = QInputDialog::getInt(this, tr("Digite o raio do corte"),
                                   tr("Limiar:"), 77, 0, dct_input.size(), 1, &ok);
    if(ok){
        QImage mod = img;
        dct_output = dct_input;

        int m = dct_input.size();
        int n = dct_input[0].size();

        for(int i=0; i < m; ++i){
            for(int j=0; j < n; ++j){
                if(i*i + j*j > cut*cut){
                    dct_output[i][j] = 0;
                    mod.setPixel(j, i, qRgb(0, 0, 0));
                }
            }
        }
        QPixmap pix = QPixmap::fromImage(mod);
        int w = ui->output_image->width();
        int h = ui->output_image->height();
        ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
        res = mod;
    }
}


void MainWindow::on_actionFiltragem_passa_alta_DCT_triggered()
{
    bool ok=false;
    int cut = QInputDialog::getInt(this, tr("Digite o raio do corte"),
                                   tr("Limiar:"), 77, 0, dct_input.size(), 1, &ok);
    if(ok){
        QImage mod = img;
        dct_output = dct_input;

        int m = dct_input.size();
        int n = dct_input[0].size();

        for(int i=0; i < m; ++i){
            for(int j=0; j < n; ++j){
                if(i*i + j*j < cut*cut){
                    dct_output[i][j] = 0;
                    mod.setPixel(j, i, qRgb(0, 0, 0));
                }
            }
        }
        QPixmap pix = QPixmap::fromImage(mod);
        int w = ui->output_image->width();
        int h = ui->output_image->height();
        ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
        res = mod;
    }
}

void MainWindow::on_actionInserir_ru_do_sal_clicando_triggered()
{
    auto mod = std::make_shared<QImage>(img); // shared copy of image
    auto dct_temp = std::make_shared<std::vector<std::vector<float>>>(dct_input); // shared DCT

    QDialog* dialog = new QDialog(this);
    dialog->setWindowTitle("Inserir sal");

    QSize displaySize(500, 500);

    ClickableLabel* label = new ClickableLabel(dialog);
    label->setFixedSize(displaySize);
    label->setAlignment(Qt::AlignCenter);
    label->setPixmap(QPixmap::fromImage(*mod).scaled(displaySize, Qt::KeepAspectRatio));
    label->setScaledContents(true);

    QPushButton* confirmButton = new QPushButton("Confirmar", dialog);

    QVBoxLayout* layout = new QVBoxLayout(dialog);
    layout->addWidget(label);
    layout->addWidget(confirmButton);
    dialog->setLayout(layout);

    connect(label, &ClickableLabel::clicked, this,
            [mod, label, displaySize, dct_temp](const QPoint& pos) mutable {
                QSize labelSize = label->size();
                QSize imgSize = mod->size();

                int x = pos.x() * imgSize.width() / labelSize.width();
                int y = pos.y() * imgSize.height() / labelSize.height();

                if (x >= 0 && x < imgSize.width() && y >= 0 && y < imgSize.height()) {
                    mod->setPixelColor(x, y, Qt::white);
                    if (!dct_temp->empty())
                        (*dct_temp)[y][x] = 255.0f;

                    QPixmap updatedPixmap = QPixmap::fromImage(*mod);
                    label->setPixmap(updatedPixmap.scaled(displaySize, Qt::KeepAspectRatio));
                }
            });

    connect(confirmButton, &QPushButton::clicked, this,
            [this, mod, dct_temp, dialog]() {
                res = *mod;
                QPixmap pix = QPixmap::fromImage(res);
                ui->output_image->setPixmap(pix.scaled(ui->output_image->size(), Qt::KeepAspectRatio));
                dct_output = *dct_temp;
                dialog->accept();
            });

    dialog->exec();
}

void MainWindow::on_actionFiltro_do_m_nimo_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int matrix[9];
            matrix[0] = qRed(img.pixel(j-1, i-1));
            matrix[1] = qRed(img.pixel(j-1, i));
            matrix[2] = qRed(img.pixel(j-1, i+1));
            matrix[3] = qRed(img.pixel(j, i-1));
            matrix[4] = qRed(img.pixel(j, i));
            matrix[5] = qRed(img.pixel(j, i+1));
            matrix[6] = qRed(img.pixel(j+1, i-1));
            matrix[7] = qRed(img.pixel(j+1, i));
            matrix[8] = qRed(img.pixel(j+1, i+1));

            std::sort(matrix, matrix + 9);

            int min = matrix[0];

            QRgb color = qRgb(min, min, min);
            mod.setPixel(j-1, i-1, color);
        }
    }
    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionFiltro_do_m_ximo_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int matrix[9];
            matrix[0] = qRed(img.pixel(j-1, i-1));
            matrix[1] = qRed(img.pixel(j-1, i));
            matrix[2] = qRed(img.pixel(j-1, i+1));
            matrix[3] = qRed(img.pixel(j, i-1));
            matrix[4] = qRed(img.pixel(j, i));
            matrix[5] = qRed(img.pixel(j, i+1));
            matrix[6] = qRed(img.pixel(j+1, i-1));
            matrix[7] = qRed(img.pixel(j+1, i));
            matrix[8] = qRed(img.pixel(j+1, i+1));

            std::sort(matrix, matrix + 9);

            int max = matrix[8];

            QRgb color = qRgb(max, max, max);
            mod.setPixel(j-1, i-1, color);
        }
    }
    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionFiltro_do_ponto_m_dio_triggered()
{
    sobelAtivo = false;
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int matrix[9];
            matrix[0] = qRed(img.pixel(j-1, i-1));
            matrix[1] = qRed(img.pixel(j-1, i));
            matrix[2] = qRed(img.pixel(j-1, i+1));
            matrix[3] = qRed(img.pixel(j, i-1));
            matrix[4] = qRed(img.pixel(j, i));
            matrix[5] = qRed(img.pixel(j, i+1));
            matrix[6] = qRed(img.pixel(j+1, i-1));
            matrix[7] = qRed(img.pixel(j+1, i));
            matrix[8] = qRed(img.pixel(j+1, i+1));

            std::sort(matrix, matrix + 9);

            int min = matrix[0];
            int max = matrix[8];

            int med = (min + max) / 2;

            QRgb color = qRgb(med, med, med);
            mod.setPixel(j-1, i-1, color);
        }
    }
    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}

QRgb grayToPseudoColor(int gray) {
    int r = 0, g = 0, b = 0;

    if (gray < 64) {
        // Preto -> Azul
        float t = gray / 64.0f;
        r = 0;
        g = 0;
        b = static_cast<int>(255 * t);
    } else if (gray < 128) {
        // Azul -> Ciano
        float t = (gray - 64) / 64.0f;
        r = 0;
        g = static_cast<int>(255 * t);
        b = 255;
    } else if (gray < 192) {
        // Ciano -> Verde
        float t = (gray - 128) / 64.0f;
        r = 0;
        g = 255;
        b = static_cast<int>(255 * (1 - t));
    } else {
        // Verde -> Amarelo
        float t = (gray - 192) / 63.0f;
        r = static_cast<int>(255 * t);
        g = 255;
        b = 0;
    }

    return qRgb(r, g, b);
}


void MainWindow::on_actionPseudo_Cores_triggered()
{
    sobelAtivo = false;
    QImage mod = img;

    int w = ui->output_image->width();
    int h = ui->output_image->height();

    for(int i = 0; i < mod.height(); ++i){
        for(int j = 0; j < mod.width(); ++j){
            QRgb color = mod.pixel(j, i);
            int gray = qGray(color);
            QRgb pseudo = grayToPseudoColor(gray);
            mod.setPixel(j, i, pseudo);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = mod;
}


void MainWindow::on_actionEqualizar_HSL_triggered()
{
    sobelAtivo = false;

    int largura = img.width();
    int altura = img.height();
    int tamanho = largura * altura;

    QVector<int> histograma(256, 0);
    QVector<QVector3D> hslPixels;

    // Converte para HSL e preenche o histograma de L
    for (int y = 0; y < altura; ++y) {
        for (int x = 0; x < largura; ++x) {
            QColor cor(img.pixel(x, y));
            float h, s, l;
            cor.getHslF(&h, &s, &l);
            int l_int = static_cast<int>(l * 255.0);
            histograma[l_int]++;
            hslPixels.append(QVector3D(h, s, l));
        }
    }

    // CDF
    QVector<int> cdf(256, 0);
    cdf[0] = histograma[0];
    for (int i = 1; i < 256; ++i)
        cdf[i] = cdf[i - 1] + histograma[i];

    QVector<int> novaL(256, 0);
    for (int i = 0; i < 256; ++i) {
        novaL[i] = qRound(((cdf[i] - cdf[0]) / double(tamanho - cdf[0])) * 255);
    }

    // Cria a nova imagem
    QImage imagemEqualizada(largura, altura, QImage::Format_RGB32);
    int idx = 0;
    for (int y = 0; y < altura; ++y) {
        for (int x = 0; x < largura; ++x) {
            QVector3D hsl = hslPixels[idx++];
            int l_int = static_cast<int>(hsl.z() * 255.0);
            float l_eq = novaL[l_int] / 255.0f;

            QColor novaCor;
            novaCor.setHslF(hsl.x(), hsl.y(), l_eq);
            imagemEqualizada.setPixel(x, y, novaCor.rgb());
        }
    }

    // Exibe na interface
    int w = ui->output_image->width();
    int h = ui->output_image->height();
    QPixmap pix = QPixmap::fromImage(imagemEqualizada);
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    res = imagemEqualizada;
}

int otsuThreshold(const QImage &grayImg) {
    int hist[256] = {0};
    int total = grayImg.width() * grayImg.height();

    // Construir histograma
    for (int y = 0; y < grayImg.height(); ++y) {
        const uchar* line = grayImg.scanLine(y);
        for (int x = 0; x < grayImg.width(); ++x) {
            int val = line[x];
            hist[val]++;
        }
    }

    // Probabilidades e média global
    float sum = 0;
    for (int t = 0; t < 256; ++t)
        sum += t * hist[t];

    float sumB = 0;
    int wB = 0, wF = 0;
    float varMax = 0;
    int threshold = 0;

    for (int t = 0; t < 256; ++t) {
        wB += hist[t];
        if (wB == 0) continue;

        wF = total - wB;
        if (wF == 0) break;

        sumB += static_cast<float>(t * hist[t]);

        float mB = sumB / wB;
        float mF = (sum - sumB) / wF;

        float varBetween = static_cast<float>(wB) * static_cast<float>(wF) * (mB - mF) * (mB - mF);

        if (varBetween > varMax) {
            varMax = varBetween;
            threshold = t;
        }
    }

    return threshold;
}

int MainWindow::otsuThreshold(const QImage &grayImg) {
    int hist[256] = {0};
    int total = grayImg.width() * grayImg.height();

    // Construir histograma
    for (int y = 0; y < grayImg.height(); ++y) {
        const uchar* line = grayImg.scanLine(y);
        for (int x = 0; x < grayImg.width(); ++x) {
            int val = line[x];
            hist[val]++;
        }
    }

    // Probabilidades e média global
    float sum = 0;
    for (int t = 0; t < 256; ++t)
        sum += t * hist[t];

    float sumB = 0;
    int wB = 0, wF = 0;
    float varMax = 0;
    int threshold = 0;

    for (int t = 0; t < 256; ++t) {
        wB += hist[t];
        if (wB == 0) continue;

        wF = total - wB;
        if (wF == 0) break;

        sumB += static_cast<float>(t * hist[t]);

        float mB = sumB / wB;
        float mF = (sum - sumB) / wF;

        float varBetween = static_cast<float>(wB) * static_cast<float>(wF) * (mB - mF) * (mB - mF);

        if (varBetween > varMax) {
            varMax = varBetween;
            threshold = t;
        }
    }

    return threshold;
}

void MainWindow::on_actionBinariza_o_OTSU_triggered()
{
    sobelAtivo = false;

    QImage gray = img;

    // Convertendo para cinza
    for (int i = 0; i < gray.height(); ++i) {
        for (int j = 0; j < gray.width(); ++j) {
            QRgb pixel = gray.pixel(j, i);
            int cinza = qGray(pixel);
            gray.setPixel(j, i, qRgb(cinza, cinza, cinza));
        }
    }

    int threshold = otsuThreshold(gray);

    QImage bin(gray.size(), QImage::Format_Grayscale8);

    for (int i = 0; i < gray.height(); ++i) {
        for (int j = 0; j < gray.width(); ++j) {
            int val = qGray(gray.pixel(j, i));
            int binVal = (val >= threshold) ? 255 : 0;
            bin.setPixel(j, i, qRgb(binVal, binVal, binVal));
        }
    }

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(QPixmap::fromImage(bin).scaled(w, h, Qt::KeepAspectRatio));
    res = bin;
}


void MainWindow::on_actionLimiariza_o_OTSU_triggered()
{
    sobelAtivo = false;

    QImage gray = img;

    // Convertendo para cinza
    for (int i = 0; i < gray.height(); ++i) {
        for (int j = 0; j < gray.width(); ++j) {
            QRgb pixel = gray.pixel(j, i);
            int cinza = qGray(pixel);
            gray.setPixel(j, i, qRgb(cinza, cinza, cinza));
        }
    }

    int threshold = otsuThreshold(gray);

    QImage lim(gray.size(), QImage::Format_Grayscale8);

    for (int i = 0; i < gray.height(); ++i) {
        for (int j = 0; j < gray.width(); ++j) {
            int val = qGray(gray.pixel(j, i));
            int result = (val >= threshold) ? val : 0;
            lim.setPixel(j, i, qRgb(result, result, result));
        }
    }

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(QPixmap::fromImage(lim).scaled(w, h, Qt::KeepAspectRatio));
    res = lim;
}

void MainWindow::on_actionComparar_laplaciano_triggered()
{
    sobelAtivo = false;
    QImage laplaciano(img.width() - 2, img.height() - 2, QImage::Format_RGB32);
    QImage laplaciano_gauss(img.width() - 4, img.height() - 4, QImage::Format_RGB32);

    // Laplaciano simples
    for(int i=1; i < img.height() - 1; ++i){
        for(int j=1; j < img.width() - 1; ++j){
            int sum = 0;
            sum += -qGray(img.pixel(j-1, i));
            sum += -qGray(img.pixel(j, i-1));
            sum +=  4*qGray(img.pixel(j, i));
            sum += -qGray(img.pixel(j, i+1));
            sum += -qGray(img.pixel(j+1, i));

            int new_color = qBound(0, sum, 255);
            laplaciano.setPixel(j-1, i-1, qRgb(new_color, new_color, new_color));
        }
    }

    // Laplaciano do Gaussiano
    for(int i=2; i < img.height() - 2; ++i){
        for(int j=2; j < img.width() - 2; ++j){
            int sum = 0;
            sum += -qGray(img.pixel(j-2, i));

            sum += -qGray(img.pixel(j-1, i-1));
            sum += -2 * qGray(img.pixel(j-1, i));
            sum += -qGray(img.pixel(j-1, i+1));

            sum += -qGray(img.pixel(j, i-2));
            sum += -2 * qGray(img.pixel(j, i-1));
            sum += 16 * qGray(img.pixel(j, i));
            sum += -2 * qGray(img.pixel(j, i+1));
            sum += -qGray(img.pixel(j, i+2));

            sum += -qGray(img.pixel(j+1, i-1));
            sum += -2 * qGray(img.pixel(j+1, i));
            sum += -qGray(img.pixel(j+1, i+1));

            sum += -qGray(img.pixel(j+2, i));

            int new_color = qBound(0, sum, 255);
            laplaciano_gauss.setPixel(j-2, i-2, qRgb(new_color, new_color, new_color));
        }
    }

    // Criação do diálogo
    QDialog* dialog = new QDialog(this);
    dialog->setWindowTitle("Comparando laplacianos");

    QSize displaySize(400, 400);

    // Imagens
    QLabel* label1 = new QLabel(dialog);
    label1->setPixmap(QPixmap::fromImage(laplaciano).scaled(displaySize, Qt::KeepAspectRatio));
    label1->setAlignment(Qt::AlignCenter);

    QLabel* label2 = new QLabel(dialog);
    label2->setPixmap(QPixmap::fromImage(laplaciano_gauss).scaled(displaySize, Qt::KeepAspectRatio));
    label2->setAlignment(Qt::AlignCenter);

    // Legendas
    QLabel* caption1 = new QLabel("Laplaciano", dialog);
    caption1->setAlignment(Qt::AlignCenter);
    QLabel* caption2 = new QLabel("Laplaciano do Gaussiano", dialog);
    caption2->setAlignment(Qt::AlignCenter);

    // Layouts individuais para imagem + legenda
    QVBoxLayout* vbox1 = new QVBoxLayout();
    vbox1->addWidget(label1);
    vbox1->addWidget(caption1);

    QVBoxLayout* vbox2 = new QVBoxLayout();
    vbox2->addWidget(label2);
    vbox2->addWidget(caption2);

    // Layout horizontal para colocar os dois conjuntos lado a lado
    QHBoxLayout* hbox = new QHBoxLayout();
    hbox->addLayout(vbox1);
    hbox->addLayout(vbox2);

    dialog->setLayout(hbox);
    dialog->exec(); // Modal
}

