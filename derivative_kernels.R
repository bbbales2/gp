# functions {
#   real QQ(real tj, real tk, real l) {
#     return exp(-((tj - tk)^2/(2 * l^2)));
#   }
#   
#   real QR(real tj, real tk, real l) {
#     return exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk) / l^2;
#   }
#   
#   real RQ(real tj, real tk, real l) {
#     return QR(tk, tj, l);
#   }
#   
#   real RR(real tj, real tk, real l) {
#     return (exp(-((tj - tk)^2/(2 * l^2))))/l^2 - (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^4;
#   }
#   
#   real QT(real tj, real tk, real l) {
#     return -((exp(-((tj - tk)^2/(2 * l^2))))/l^2) + (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^4;
#   }
#   
#   real TQ(real tj, real tk, real l) {
#     return QT(tk, tj, l);
#   }
#   
#   real RT(real tj, real tk, real l) {
#     return (3 * exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk))/l^4 - (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^3)/l^6;
#   }
#   
#   real TR(real tj, real tk, real l) {
#     return RT(tk, tj, l);
#   }
#   
#   real TT(real tj, real tk, real l) {
#     return (3 * exp(-((tj - tk)^2/(2 * l^2))))/l^4 - (6 * exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^6 + (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^4)/l^8;
#   }
# }

QQ = function(tj, tk, l) {
  exp(-((tj - tk)^2/(2 * l^2)))
}

QR = function(tj, tk, l) {
  (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk))/l^2
}

RQ = function(tj, tk, l) {
  QR(tk, tj, l)
}

RR = function(tj, tk, l) {
  (exp(-((tj - tk)^2/(2 * l^2))))/l^2 - (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^4
}

QT = function(tj, tk, l) {
  -((exp(-((tj - tk)^2/(2 * l^2))))/l^2) + (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^4
}

TQ = function(tj, tk, l) {
  QT(tk, tj, l)
}

RT = function(tj, tk, l) {
  (3 * exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk))/l^4 - (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^3)/l^6
}

TR = function(tj, tk, l) {
  RT(tk, tj, l)
}

TT = function(tj, tk, l) {
  (3 * exp(-((tj - tk)^2/(2 * l^2))))/l^4 - (6 * exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^2)/l^6 + (exp(-((tj - tk)^2/(2 * l^2))) * (tj - tk)^4)/l^8
}
