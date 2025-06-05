#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QImage>
#include <QPixmap>
#include <QRandomGenerator>
#include <QDebug>
#include <QInputDialog>
#include <QtMath>
#include <QtConcurrent/QtConcurrent>
#include "clickablelabel.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    QImage sobelImage;
    std::vector<std::vector<int>> magnitudes;
    std::vector<std::vector<float>> angles;
    bool sobelAtivo = false;

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;

private slots:
    void on_actionEscala_Cinza_triggered();

    void on_actionCarregar_imagem_triggered();

    void on_actionEscala_Cinza_Inversa_triggered();

    void on_actionColorida_Inversa_triggered();

    void on_actionSeparar_canais_de_cores_triggered();

    void on_actionSal_e_pimenta_triggered();

    void on_output_to_input_btn_clicked();

    void on_actionFiltro_mediana_triggered();

    void on_actionFiltro_m_dia_triggered();

    void on_actionBinarizar_triggered();

    void on_actionLaplaciano_triggered();

    void on_actionEqualiza_o_triggered();

    void on_actionCanal_Vermelho_em_Cinza_triggered();

    void on_actionRGB_para_HSV_triggered();

    void on_actionHSV_para_RGB_triggered();

    void on_actionCompress_o_de_Escala_Din_mica_triggered();

    void on_actionBordas_por_Sobel_triggered();

    void on_actionLimiariza_o_triggered();

    void on_actionDCT_triggered();

    void on_actionIDCT_triggered();

    void on_actionFiltragem_passa_baixa_DCT_triggered();

    void on_actionFiltragem_passa_alta_DCT_triggered();

    void on_actionInserir_ru_do_sal_clicando_triggered();

    void on_actionFiltro_do_m_nimo_triggered();

    void on_actionFiltro_do_m_ximo_triggered();

    void on_actionFiltro_do_ponto_m_dio_triggered();

    void on_actionPseudo_Cores_triggered();

    void on_actionEqualizar_HSL_triggered();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
