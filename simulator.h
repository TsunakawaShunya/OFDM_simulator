/*
 * File:   simulator.h
 * Author: Tsunakawa
 *
 * Created on 2024/12/20, 18:31
*/

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES
#include </Volumes/USB1/eigen-3.4.0/Eigen/Core>
#include </Volumes/USB1/eigen-3.4.0/Eigen/Eigen>
#include <cmath>
#include <utility>
#include "random_collection.h"

class Simulator {
    public:
        Simulator(double weightOfLastPath) {
            // 試行回数を設定
            std::cout << "--------------------------------------------------------------------" << std::endl;
            std::cout << "TRIAL?" << std::endl;
            std::cout << "--------------------------------------------------------------------" << std::endl;
            std::cin >> NUMBER_OF_TRIAL;

            // リサイズ
            W_.resize(Q_, K_);
            quzaiSqrt_.resize(Q_);
            h_.resize(L_, Q_);
            H_.resize(L_, K_);
            H_est_.resize(K_);
            txData_.resize(L_, K_);
            rxData_.resize(L_, K_);
            Y_.resize(L_, K_);
            X_.resize(L_, K_);
            R_.resize(L_, K_);
            Cmat_.resize(L_, L_);
            symbol_.resize(NUMBER_OF_SYMBOLS);

            // DFT行列設定
            setW_();

            // 共分散行列設定
            setCmat_();
            // シンボル設計
            setSymbol();

            // 乱数設定
            unitIntUniformRand_.init(0, NUMBER_OF_SYMBOLS - 1, seed);
            unitCNormalRand_.init(0.0, 1 / sqrt(2), seed);

            // 遅延プロファイル生成
            setChannelProfile(weightOfLastPath);
        }

        virtual ~Simulator() {
        }


        /**
         * 雑音の分散設定
         * @param EbN0dB EbN0 [dB]
         */
        void setNoiseSD(double EbN0dB) {
            noiseSD_ = std::sqrt(std::pow(10.0, -0.1 * EbN0dB) / (double)NUMBER_OF_BIT);
        }

        /**
         * ドップラー周波数設定
         * @param f_d ドップラー周波数
         */
        void setDopplerFrequence(double f_d) {
            f_d_ = f_d;
        }

        /**
         * 数値計算実験
         * @return ビット誤り率のシミュレーション値
         */
        double getBERSimulation() {
            int count = 0;
            for(auto tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setX_();
                setH_();
                setY_();
                estimateChannelByML();
                equalizeByEstimatedChannel();
                setRxDataByML();
                count += getBitErrorCount();
            }
            return (double)count / (double)NUMBER_OF_TRIAL / (double)NUMBER_OF_BIT / (double)K_ / (double)L_;
        }

        /**
         * 数値計算実験
         * @return 平均二乗誤差のシミュレーション値
         */
        double getMSESimulation() {
            double mse = 0.0;
            for(auto tri = 0; tri < NUMBER_OF_TRIAL; tri++) {
                setX_();
                setH_();
                setY_();
                estimateChannelByML();
                mse += getMeanSquaredError();
            }
            return mse / NUMBER_OF_TRIAL;
        }

    private:
        double noiseSD_;                            // 雑音の標準偏差
        int NUMBER_OF_TRIAL;                        // 試行回数
        const int K_ = 52;                          // サブキャリア数 K
        const int L_ = 10;                          // 1フレームのシンボル数（わからん） L
        const int Q_ = 16;                          // 伝送路のインパルス応答のパス数 Q
        const double T_ = 32.0;                     // 有効シンボル長 T
        const double TGI_ = 8.0;                    // ガードインターバル長 TGI
        double Ts_ = T_ + TGI_;                     // シンボル長
        const int NUMBER_OF_FFT = K_;               // FFTポイント数(IEEE802.11aとかは64だった気がする)
        const int NUMBER_OF_PILOT = 2;              // パイロットシンボル個数
        const int NUMBER_OF_SYMBOLS = 2;            // 変調方式に合わせて(BPSKの例)
        const int NUMBER_OF_BIT = 1;                // 変調方式に合わせて(BPSKの例)

        double decayConstant_;                      // 伝送路のインパルス応答の指数関数モデルの減衰定数
        double f_d_;                                // ドップラー周波数

        Eigen::MatrixXcd W_;            // DFT行列:式(17)
        Eigen::VectorXd quzaiSqrt_;     // 伝送路のインパルス応答の遅延プロファイル（ルート）
        Eigen::MatrixXcd h_;            // インパルス応答
        Eigen::MatrixXcd H_;            // 周波数応答
        Eigen::MatrixXi txData_;        // 送信アルファベット
        Eigen::MatrixXi rxData_;        // 受信アルファベット
        Eigen::VectorXcd symbol_;       // 送信可能なシンボルベクトル
        Eigen::MatrixXcd Y_;            // 受信信号
        Eigen::MatrixXcd X_;            // 送信信号
        Eigen::MatrixXcd R_;            // 等化後の受信信号
        Eigen::MatrixXcd Cmat_;         // Jakesモデルの共分散行列
        Eigen::VectorXcd H_est_;        // 伝送路の周波数応答の推定値

        // 乱数用
        int seed = 100;
        uniform_int_distribution<> unitIntUniformRand_;     // int型一様乱数
        cnormal_distribution<> unitCNormalRand_;            // 平均0，分散1の複素正規分布（実部，虚部それぞれ平均0，分散0.5の正規分布）

        // DFT行列Wの生成:式(17)
        void setW_() {
            for (auto q = 0; q < Q_; ++q) {
                for (auto k = 0; k < K_ / 2; ++k) {
                    // -26から-1番目のキャリヤ
                    W_(q, k) = std::polar(1.0, -2.0 * M_PI * (k - K_ / 2) * q / NUMBER_OF_FFT);
                    // 1から26番目のキャリヤ
                    W_(q, k + K_ / 2) = std::polar(1.0, -2.0 * M_PI * (k + 1) * q / NUMBER_OF_FFT);
                }
            }
        }

        /**
         * 遅延プロファイル生成
         * @param 最後のパスの重み
         */
        void setChannelProfile(double weightOfLastPath) {
            assert(0 < weightOfLastPath && weightOfLastPath <= 1.0);
            // 伝送路プロファイルの減衰係数
            decayConstant_ = - NUMBER_OF_FFT * std::log(weightOfLastPath) / T_ / (Q_ - 1);

            // 伝送路プロファイルの生成
            for (auto q = 0; q < Q_; ++q) {
                quzaiSqrt_(q) = std::exp(-decayConstant_ * q);
            }
            auto tmp = quzaiSqrt_.sum();
            // プロファイルの和を1に正規化してルート:式(20)
            for (auto q = 0; q < Q_; ++q) {
                quzaiSqrt_(q) = std::sqrt(quzaiSqrt_(q) / tmp);
            }
        }

        /**
         * 送信信号生成
         */
        void setX_() {
            for(auto k = 0; k < K_; k++) {
                for(auto l = 0; l < L_; l++) {
                    if(l < NUMBER_OF_PILOT) {
                        txData_(l, k) = 1;      // パイロットシンボルのデータわからない
                    } else {
                        txData_(l, k) = unitIntUniformRand_();
                    }
                    X_(l, k) = symbol_(txData_(l, k));
                }
            }
        }

        /**
         * 周波数応答生成
         */
        void setH_() {
            // 伝送路のインパルス応答の生成
            // おそらくここに時間変化のJakesモデルを入れる
            seth_();
            // 伝送路の周波数応答の生成:式(19)
            H_ = h_ * W_;
        }

        /**
         * Jakesモデルの共分散行列生成
         */
        void setCmat_() {
            for(auto l_1 = 0; l_1 < L_; l_1++) {
                for(auto l_2 = 0; l_2 < L_; l_2++) {
                    Cmat_(l_1, l_2) = j0(2 * M_PI * abs((l_1 - l_2) * Ts_) * f_d_);
                }
            }
        }

        /**
         * インパルス応答行列生成
         */
        void seth_() {
            Eigen::VectorXcd h_q(L_);
            Eigen::MatrixXcd A(L_, L_);
            Eigen::VectorXcd x(L_);

            // Cmat_ = U * Lambda * U^-1 に固有値分解してUとΛを取得
            auto [U, Lambda] = EigenvalueDecomposition(Cmat_);
            // A生成
            A = U * Lambda.array().sqrt().matrix();


            for(auto q = 0; q < Q_; q++) {
                for(auto l = 0; l < L_; l++) {
                    x(l) = quzaiSqrt_(q) * unitCNormalRand_();
                    h_q = A * x;        // h_q = Ax
                    h_(l, q) = h_q(l);
                }
            }
        }

        /**
         * 固有値分解
         * @param mat 固有値分解する行列
         * @return std::pair ユニタリー行列U, 固有値行列Lambda
         */
        std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> EigenvalueDecomposition(const Eigen::MatrixXcd& mat) {
            Eigen::EigenSolver<Eigen::MatrixXcd> solver(mat);

            // 固有値を取得 (対角行列として表現)
            Eigen::MatrixXcd Lambda = Eigen::MatrixXcd::Zero(L_, L_);
            Eigen::VectorXcd eigenvalues = solver.eigenvalues(); // 固有値ベクトル
            for (int i = 0; i < L_; ++i) {
                Lambda(i, i) = eigenvalues(i); // 固有値を対角行列に配置
            }

            // 固有ベクトル (ユニタリー行列 U) を取得
            Eigen::MatrixXcd U = solver.eigenvectors();

            return std::make_pair(U, Lambda);
        }

        /**
         * 受信信号生成
         */
        void setY_() {
            for (auto l = 0; l < L_; l++) {
                for (auto k = 0; k < K_; k++) {
                    Y_(l, k) = H_(l, k) * X_(l, k) + noiseSD_ * unitCNormalRand_();
                }
            }
        }

        /**
         * 平均二乗誤差（L2ノルム）
         * @return 誤差のL2ノルム
         */
        double getMeanSquaredError() {
            double sum = 0;
            for(auto k = 0; k < K_; k++) {
                for(auto l = 0; l < L_; l++) {
                    sum += std::norm(H_(l, k) - H_est_(k));
                }
            }
            return sqrt(sum / (double)L_ / (double)K_);
        }

        /**
         * MLによる伝送路推定
         */
        void estimateChannelByML() {
            Eigen::VectorXcd sum(K_);
            sum.setZero();
            for(auto k = 0; k < K_; k++) {
                for(auto i = 0; i < NUMBER_OF_PILOT; i++) {
                    // パイロットシンボル区間での推定
                    sum(k) += Y_(k, i) / X_(k, i);
                }
            }
            H_est_ = sum / NUMBER_OF_PILOT;
        }

        /**
         * 推定値を用いた等化
         */
        void equalizeByEstimatedChannel() {
            for(auto l = 0; l < L_; l++) {
                for(auto k = 0; k < K_; k++) {
                    R_(l, k) = Y_(l, k) / H_est_(k);
                }
            }
        }

        /**
         * 完璧な等化
         */
        void equalizePerfectly() {
            for(auto l = 0; l < L_; l++) {
                for(auto k = 0; k < K_; k++) {
                    R_(l, k) = Y_(l, k) / H_(l, k);
                }
            }
        }


        /**
         * シンボル生成（BPSK）
         */
        void setSymbol() {
            symbol_(0) = -1;
            symbol_(1) = 1;
        }

        /**
         * 最尤復調
         */
        void setRxDataByML() {
            Eigen::VectorXd obj(NUMBER_OF_SYMBOLS);     // 最小化の目的関数

            for(auto l = 0; l < L_; l++) {
                for(auto k = 0; k < K_; k++) {
                    for(auto i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                        // 最尤復調の周波数応答は推定値を使う？
                        obj(i) = std::norm(H_est_(k)) * std::norm((R_(l, k) - symbol_(i)));
                    }
                    Eigen::VectorXd::Index minColumn;       // ノルムが最小な index（つまり受信データ）
                    obj.minCoeff(&minColumn);
                    rxData_(l, k) = minColumn;
                }
            }
        }

        /**
         * 最尤復調（完璧な推定値の場合）
         */
        void setRxDataByPerfect() {
            Eigen::VectorXd obj(NUMBER_OF_SYMBOLS);     // 最小化の目的関数

            for(auto l = 0; l < L_; l++) {
                for(auto k = 0; k < K_; k++) {
                    for(auto i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                        obj(i) = std::norm(H_(l, k)) * std::norm((R_(l, k) - symbol_(i)));
                    }
                    Eigen::VectorXd::Index minColumn;       // ノルムが最小な index（つまり受信データ）
                    obj.minCoeff(&minColumn);
                    rxData_(l, k) = minColumn;
                }
            }
        }

        /**
         * ビット誤り数のカウント
         * @return 全ての誤りビット数
         */
        int getBitErrorCount() {
            int count = 0;
            for(auto l = 0; l < L_; l++) {
                for(auto k = 0; k < K_; k++) {
                    count += hammingDistance(txData_(l, k), rxData_(l, k));
                }
            }
            return count;
        }

        /**
         * ハミング距離計算
         * @param 整数1，整数2
         * @return ハミング距離
         */
        int hammingDistance(int num1, int num2) {
            int ham = 0;
            int xorResult;
            int bitMask = 1;

            xorResult = num1 ^ num2;

            for(int i = 0; i < NUMBER_OF_BIT; i++) {
                ham += (xorResult & bitMask) >> i;
                bitMask <<= 1;
            }
            
            return ham;
        }
};

#endif /* SIMULATOR_H */