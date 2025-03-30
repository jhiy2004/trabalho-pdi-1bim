#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QImage>
#include <QPixmap>
#include <QRandomGenerator>
#include <QDebug>

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

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
