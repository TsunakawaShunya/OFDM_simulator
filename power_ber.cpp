/*
 * File:   power_ber.cpp
 * Author: Tsunakawa
 *
 * Created on 2024/12/20, 19:36
 */
#include <iostream>
#include <fstream>
#include "simulator.h"

// パラメータ
static const double WEIGHT_OF_LAST_PATH = 0.1;
static const int EbN0dBmin = 0;
static const int EbN0dBmax = 30;
static const int EbN0dBstp = 5;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  

// BER
double ber;

// ドップラー周波数
double dopplerFrequence;

int main()
{
    Simulator sim(WEIGHT_OF_LAST_PATH);

    // ドップラー周波数を設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "f_d?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> dopplerFrequence;
    sim.setDopplerFrequence(dopplerFrequence);

    fileName = "BER_f_d=" + std::to_string(dopplerFrequence) + ".csv";
	ofs.open(fileName);


    for (auto EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        // SN設定
        sim.setNoiseSD(EbN0dB);
        // シミュレーション
        ber = sim.getBERSimulation();
        // 標準出力
        std::cout << EbN0dB << "," << ber << std::endl;
        // CSV出力
		ofs << EbN0dB << "," << ber << std::endl;
    }
    ofs.close();
    
    return 0;
}