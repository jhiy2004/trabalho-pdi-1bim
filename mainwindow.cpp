#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QMessageBox>
#include <QMouseEvent>

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

void MainWindow::on_actionEscala_Cinza_triggered()
{
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
}

void MainWindow::on_actionCarregar_imagem_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open a file", ".");
    img.load(filename);

    int w = ui->input_image->width();
    int h = ui->input_image->height();
    QPixmap pix(filename);
    ui->input_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
}

void MainWindow::on_actionEscala_Cinza_Inversa_triggered()
{
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
}


void MainWindow::on_actionColorida_Inversa_triggered()
{
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
}


void MainWindow::on_actionSeparar_canais_de_cores_triggered()
{
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
}

void MainWindow::on_output_to_input_btn_clicked()
{
    QPixmap output_pix = ui->output_image->pixmap();
    if (!output_pix || output_pix.isNull())
        return;

    img = output_pix.toImage();
    QPixmap scaledPix = output_pix.scaled(output_pix.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation);
    ui->input_image->setPixmap(scaledPix);
}



void MainWindow::on_actionFiltro_mediana_triggered()
{
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
}

void MainWindow::on_actionFiltro_m_dia_triggered()
{
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
}


void MainWindow::on_actionBinarizar_triggered()
{
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
    }
}

void MainWindow::on_actionLaplaciano_triggered()
{
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    // Garante que o formato da imagem é compatível com setPixel()
    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    std::vector<std::vector<int>> magnitudes;
    int min_magnitude = INT_MAX;
    int max_magnitude = INT_MIN;

    for(int i=1; i < img.height() - 1; ++i){
        std::vector<int> row;
        for(int j=1; j < img.width() - 1; ++j){
            int sum=0;
            sum += -qGray(img.pixel(j-1, i));
            sum += -qGray(img.pixel(j, i-1));
            sum += 4*qGray(img.pixel(j, i));
            sum += -qGray(img.pixel(j, i+1));
            sum += -qGray(img.pixel(j+1, i));

            row.push_back(sum);

            if (sum < min_magnitude) min_magnitude = sum;
            if (sum > max_magnitude) max_magnitude = sum;
        }
        magnitudes.push_back(row);
    }

    for(int i=0; i < img.height() - 2; ++i){
        for(int j=0; j < img.width() - 2; ++j){
            int magnitude = magnitudes[i][j];
            int norm = 0;
            if (max_magnitude != min_magnitude) {
                norm = (magnitude - min_magnitude) * 255 / (max_magnitude - min_magnitude);
            }
            QRgb color = qRgb(norm, norm, norm);
            mod.setPixel(j, i, color);
        }
    }

    QPixmap pix = QPixmap::fromImage(mod);

    int w = ui->output_image->width();
    int h = ui->output_image->height();
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
}


void MainWindow::on_actionEqualiza_o_triggered()
{
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
}


void MainWindow::on_actionCanal_Vermelho_em_Cinza_triggered()
{
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
        H += 360;
    }

    // Transformo em porcentagem
    S = S * 100;
    V = V * 100;

    QString result = QString("HSV:\nH: %1\nS: %2\nV: %3").arg(H).arg(S).arg(V);
    QMessageBox::information(this, "Resultado HSV", result);
}


void MainWindow::on_actionHSV_para_RGB_triggered()
{
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
}


void MainWindow::on_actionBordas_por_Sobel_triggered()
{
    QImage mod(img.width() - 2, img.height() - 2, QImage::Format_RGB32);

    if (mod.format() != QImage::Format_RGB32 && mod.format() != QImage::Format_ARGB32) {
        mod = mod.convertToFormat(QImage::Format_RGB32);
    }

    int min_magnitude = INT_MAX;
    int max_magnitude = INT_MIN;

    for(int i = 1; i < img.height() - 1; ++i){
        std::vector<int> row;
        for(int j = 1; j < img.width() - 1; ++j){
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
            row.push_back(magnitude);

            if (magnitude < min_magnitude) min_magnitude = magnitude;
            if (magnitude > max_magnitude) max_magnitude = magnitude;
        }
        magnitudes.push_back(row);
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
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (obj == ui->output_image && event->type() == QEvent::MouseMove && sobelAtivo) {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);

        if (sobelImage.isNull()) return true;

        QPoint pos = mouseEvent->pos();

        // Tamanho do label e da imagem
        QSize labelSize = ui->output_image->size();
        QSize imageSize = sobelImage.size();

        // Calcula escala mantendo proporção
        double scale = qMin(
            static_cast<double>(labelSize.width()) / imageSize.width(),
            static_cast<double>(labelSize.height()) / imageSize.height()
            );

        // Tamanho da imagem escalada
        QSize scaledSize = imageSize * scale;

        // Offset da imagem dentro do QLabel (centralizado)
        int offsetX = (labelSize.width() - scaledSize.width()) / 2;
        int offsetY = (labelSize.height() - scaledSize.height()) / 2;

        // Verifica se o mouse está dentro da área da imagem desenhada
        if (pos.x() < offsetX || pos.x() >= offsetX + scaledSize.width() ||
            pos.y() < offsetY || pos.y() >= offsetY + scaledSize.height()) {
            ui->magnitude_label->clear();
            return true;
        }

        // Mapeia posição do mouse para coordenada na imagem original
        int imgX = (pos.x() - offsetX) * imageSize.width() / scaledSize.width();
        int imgY = (pos.y() - offsetY) * imageSize.height() / scaledSize.height();

        // Verifica limites e exibe magnitude
        if (imgX >= 0 && imgX < imageSize.width() &&
            imgY >= 0 && imgY < imageSize.height()) {

            int mag = qGray(sobelImage.pixel(imgX, imgY));
            ui->magnitude_label->setText(QString("Magnitude: %1").arg(mag));
        } else {
            ui->magnitude_label->clear();
        }

        return true;
    }

    return QMainWindow::eventFilter(obj, event);
}


void MainWindow::on_actionLimiariza_o_triggered()
{
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
    }
}
