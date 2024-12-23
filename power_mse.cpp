
/*
 * File:   power_mse.cpp
 * Author: Tsunakawa
 *
 * Created on 2024/12/20, 19:36
 */
#include <iostream>
#include <fstream>
#include "simulator.h"

// パラメータ
static const double WEIGHT_OF_LAST_PATH = 0.1;
static const int dopplerFrequence_min = 0;
static const int dopplerFrequence_max = 1;
static const int dopplerFrequence_stp = 0.1;

// ファイル
std::string fileName;       // ファイル名
std::ofstream ofs;        // 出力ファイル  

// MSE
double mse;

// SN
double EbN0dB;

int main() {
	fileName = "MSE.csv";
	ofs.open(fileName);

    Simulator sim(WEIGHT_OF_LAST_PATH);

    // SNを設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Eb/N0?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> EbN0dB;
    sim.setNoiseSD(EbN0dB);

    for (auto dopplerFrequence = dopplerFrequence_min; dopplerFrequence <= dopplerFrequence_max; dopplerFrequence += dopplerFrequence_stp) {
        // SN設定
        sim.setNoiseSD(EbN0dB);
        // シミュレーション
        mse = sim.getMSESimulation();
        // 標準出力
        std::cout << "-----------" << std::endl;
        std::cout << EbN0dB << "," << mse << std::endl;
        // CSV出力
		ofs << EbN0dB << "," << mse << std::endl;
    }
    ofs.close();

    return 0;
}