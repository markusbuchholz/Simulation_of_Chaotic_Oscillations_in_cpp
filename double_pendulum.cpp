#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//----------- system dynamic parameters --------------------

float m1 = 2.0;
float m2 = 2.5;
float l1 = 1.5; //23.0;
float l2 = 2.0; //16.0;
float g = 9.81;


float dt = 0.001;

//-----------------------------------------------------------
// theta1_dot
float function1(float theta1, float theta2, float theta1_dot, float theta2_dot)
{

    return theta1_dot;
}

//-----------------------------------------------------------
// theta1_dot_dot
float function2(float theta1, float theta2, float theta1_dot, float theta2_dot)
{
    //ok
    //float theta1_dot_dot = (omega2 * l * (-std::sin(theta1) + M * cos(theta1 - theta2) * sin(theta2)) - M * l * (theta1_dot * theta1_dot * cos(theta1 - theta2) + l * theta2_dot * theta2_dot) * sin(theta1 - theta2)) / (l - M * l * (std::cos(theta1 - theta2)) * (std::cos(theta1 - theta2)));
    
    
    float theta1_dot_dot = -g * (2 * m1 + m2) * std::sin(theta1) - m2 * g * std::sin(theta1 - 2 * theta2) - 2 * std::sin(theta1 - theta2) * m2 * (l2 * theta2_dot * theta2_dot + l1 * theta2_dot * theta2_dot * std::cos(theta1 - theta2)) / (l1 * (2 * m1 + m2 - m2 * std::cos(2 * theta1 - 2 * theta2)));
    return theta1_dot_dot;
}

//-----------------------------------------------------------
// theta2_dot
float function3(float theta1, float theta2, float theta1_dot, float theta2_dot)
{

    return theta2_dot;
}

//-----------------------------------------------------------
// theta2_dot_dot
float function4(float theta1, float theta2, float theta1_dot, float theta2_dot)
{
    //ok
    //float theta2_dot_dot = (omega2 * std::cos(theta1 - theta2) * std::sin(theta1) - omega2 * sin(theta2) + (theta1_dot * theta1_dot + M * l * theta2_dot * theta2_dot * std::cos(theta1 - theta2)) * std::sin(theta1 - theta2)) / (l - M * l * (std::cos(theta1 - theta2)) * (std::cos(theta1 - theta2)));
   
    float theta2_dot_dot = 2 * std::sin(theta1 - theta2) * (l1 * theta1_dot * theta1_dot * (m1 + m2) + g * (m1 + m2) * std::cos(theta1) + l2 * theta2_dot * theta2_dot * m2 * std::cos(theta1 - theta2)) / (l2 * (2 * m1 + m2 - m2 * std::cos(2 * theta1 - 2 * theta2)));

    return theta2_dot_dot;
}

//-----------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> methodRuneKuttaDoublePndulum()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> diffEq3;
    std::vector<float> diffEq4;

    std::vector<float> time;

    // init values
    float x1 = (M_PI / 180) * 40; // theta1 - pendulum
    float x2 = (M_PI / 180) * 40; // theta2 - pednulum
    float x3 = 0;                 // theta1_dot
    float x4 = 0;                 // theta2_dot
    float t = 0.0;                // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    diffEq3.push_back(x3);
    diffEq4.push_back(x4);
    time.push_back(t);

    for (int ii = 0; ii < 10000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2, x3, x4);
        float k12 = function2(x1, x2, x3, x4);
        float k13 = function3(x1, x2, x3, x4);
        float k14 = function4(x1, x2, x3, x4);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k23 = function3(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k24 = function4(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k33 = function3(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k34 = function4(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k43 = function3(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k44 = function4(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        x3 = x3 + dt / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        x4 = x4 + dt / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        diffEq3.push_back(x3);
        diffEq4.push_back(x4);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, diffEq3, diffEq4, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::tuple<std::vector<float>, std::vector<float>> data1)
{

    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> yY1 = std::get<1>(data1);

    plt::plot(xX1, yY1);
    plt::xlabel("t[s]");
    plt::ylabel("[rad]");
    plt::show();
}

//---------------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>> computeXYpos1(std::vector<float> angle)
{

    std::vector<float> xX;
    std::vector<float> yY;

    for (auto &ii : angle)
    {

        xX.push_back(std::sin(ii) * l1);
        yY.push_back(std::cos(ii) * l1);
    }

    return std::make_tuple(xX, yY);
}
//---------------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>> computeXYpos2(std::vector<float> angle1, std::vector<float> angle2)
{

    std::vector<float> xX;
    std::vector<float> yY;

    for (int ii = 0; ii < angle1.size(); ii++)
    {

        xX.push_back(std::sin(angle1[ii]) * l1 + std::sin(angle2[ii]) * l2);
        yY.push_back(std::cos(angle1[ii]) * l1 + std::cos(angle2[ii]) * l2);
    }

    return std::make_tuple(xX, yY);
}
//---------------------------------------------------------------

void plot2D2D(std::tuple<std::vector<float>, std::vector<float>> data1, std::tuple<std::vector<float>, std::vector<float>> data2)
{

    std::vector<float> xX1 = std::get<0>(data1);
    std::vector<float> yY1 = std::get<1>(data1);

    std::vector<float> xX2 = std::get<0>(data2);
    std::vector<float> yY2 = std::get<1>(data2);

    plt::plot(xX1, yY1);
    plt::plot(xX2, yY2);
    plt::xlabel("X[m]");
    plt::ylabel("Y[m]");
    plt::show();
}

//---------------------------------------------------------------
int main()
{

    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> dyn = methodRuneKuttaDoublePndulum();

    std::tuple<std::vector<float>, std::vector<float>> theta1 = std::make_tuple(std::get<4>(dyn), std::get<0>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> theta2 = std::make_tuple(std::get<4>(dyn), std::get<1>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> theta12 = std::make_tuple(std::get<0>(dyn), std::get<1>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> theta1_dot = std::make_tuple(std::get<4>(dyn), std::get<2>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> theta2_dot = std::make_tuple(std::get<4>(dyn), std::get<3>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> XYpendulum1 = computeXYpos1(std::get<0>(dyn));
    std::tuple<std::vector<float>, std::vector<float>> XYpendulum2 = computeXYpos2(std::get<0>(dyn), std::get<1>(dyn));
   //  plot2D2D(theta1, theta2);
     plot2D2D(XYpendulum1, XYpendulum2);
    // plot2D(XYpendulum1);
    // plot2D(XYpendulum2);
    plot2D(theta2);
}
