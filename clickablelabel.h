#ifndef CLICKABLELABEL_H
#define CLICKABLELABEL_H

#include <QLabel>
#include <QMouseEvent>

class ClickableLabel : public QLabel {
    Q_OBJECT
public:
    explicit ClickableLabel(QWidget* parent = nullptr) : QLabel(parent) {}
    explicit ClickableLabel(const QString& text, QWidget* parent = nullptr) : QLabel(text, parent) {}

signals:
    void clicked(const QPoint& pos);

protected:
    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            emit clicked(event->pos());
        }
        QLabel::mousePressEvent(event);
    }
};

#endif // CLICKABLELABEL_H
