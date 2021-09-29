#pragma once

#include "Eigen.h"

void print_mtxf(const Eigen::MatrixXf& X)
{
    int i, j, nrow, ncol;

    nrow = X.rows();
    ncol = X.cols();

    for (i=0; i<nrow; i++)
    {
        for (j=0; j<ncol; j++)
        {
            Serial.print(X(i,j), 3);   // print 6 decimal places
            Serial.print(" ");
        }
        Serial.println();
    }
    Serial.println();
}

void print_visual(const Eigen::MatrixXf& quat) {
  Serial.print('w');
  Serial.print(quat(0));
  Serial.print("wa");
  Serial.print(quat(1));
  Serial.print("ab");
  Serial.print(quat(2));
  Serial.print("bc");
  Serial.print(quat(3));
  Serial.println("c");
}
