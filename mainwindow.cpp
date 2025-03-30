#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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
        img.setPixel(x, y, new_c);
    }

    QPixmap pix = QPixmap::fromImage(img);
    ui->output_image->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
}

void MainWindow::on_output_to_input_btn_clicked()
{
    QPixmap output_pix = ui->output_image->pixmap();
    int w = output_pix.width();
    int h = output_pix.height();

    img = output_pix.toImage();
    ui->input_image->setPixmap(output_pix.scaled(w, h, Qt::KeepAspectRatio));
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

