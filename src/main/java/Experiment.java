import Jama.Matrix;

import java.util.Random;

public class Experiment {
    float x_min1, x_min2, x_max1, x_max2, x_min3, x_max3;
    float y_min;
    float y_max;
    double Sb_2;
    int n = 15;
    int m = 3;
    private static final int constanta = 5000;
    double[][] matrix_x = new double[n][10];
    MatrixFunc matrixFunc = new MatrixFunc();


    public Experiment(float x_min1, float x_min2, float x_min3, float x_max1, float x_max2, float x_max3) {
        this.x_min1 = x_min1;
        this.x_min2 = x_min2;
        this.x_min3 = x_min3;
        this.x_max1 = x_max1;
        this.x_max2 = x_max2;
        this.x_max3 = x_max3;
        this.y_min = 200 + (x_min1 + x_min2 + x_min3) / 3;
        this.y_max = 200 + (x_max1 + x_max2 + x_max3) / 3;

    }

    private float randFloat(float min, float max) {
        Random random = new Random();
        return random.nextFloat() * (max - min) + min;
    }

    public double[][] genMatrix() {
        double[][] matrix = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                matrix[i][j] = randFloat(y_min, y_max);
            }
        }
        return matrix;
    }


    double[][] mainMatrix = {
            {-1, -1, -1, 1, 1, 1, -1, 1, 1, 1},
            {-1, -1, 1, 1, -1, -1, 1, 1, 1, 1},
            {-1, 1, -1, -1, 1, -1, 1, 1, 1, 1},
            {-1, 1, 1, -1, -1, 1, -1, 1, 1, 1},
            {1, -1, -1, -1, -1, 1, 1, 1, 1, 1},
            {1, -1, 1, -1, 1, -1, -1, 1, 1, 1},
            {1, 1, -1, 1, -1, -1, -1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {-1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0},
            {1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0},
            {0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0},
            {0, 1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0},
            {0, 0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623},
            {0, 0, 1.215, 0, 0, 0, 0, 0, 0, 1.4623},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    };

    double[][] CohrenTable = {
            {.0, .0, .0, .0, .0, .0, .0, .0, .0, .0},
            {.9985, .9750, .9392, .9057, .8772, .8534, .8332, .8159, .8010, .7880},
            {.9669, .8709, .7977, .7457, .7071, .6771, .6530, .6333, .6167, .6025},
            {.9065, .7679, .6841, .6287, .5892, .5598, .5365, .5175, .5017, .4884},
            {.8412, .6838, .5981, .5440, .5063, .4783, .4564, .4387, .4241, .4118},
            {.7808, .6161, .5321, .4803, .4447, .4184, .3980, .3817, .3682, .3568},
            {.7271, .5612, .4800, .4307, .3974, .3726, .3535, .3384, .3259, .3154},
            {.6798, .5157, .4377, .3910, .3595, .3362, .3185, .3043, .2926, .2829},
            {.6385, .4775, .4027, .3584, .3286, .3067, .2901, .2768, .2659, .2568},
            {.6020, .4450, .3733, .3311, .3029, .2823, .2666, .2541, .2439, .2353},
            {.5410, .3924, .3264, .2880, .2624, .2439, .2299, .2187, .2098, .2020},
            {.5410, .3924, .3264, .2880, .2624, .2439, .2299, .2187, .2098, .2020},
            {.5410, .3924, .3264, .2880, .2624, .2439, .2299, .2187, .2098, .2020},
            {.4709, .3346, .2758, .2419, .2159, .2034, .1911, .1815, .1736, .1671},
            {.4709, .3346, .2758, .2419, .2159, .2034, .1911, .1815, .1736, .1671},
            {.4709, .3346, .2758, .2419, .2159, .2034, .1911, .1815, .1736, .1671},
    };

    double[] StudentaTable = {12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262,
            2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.12, 2.11, 2.101, 2.093, 2.086,
            2.08, 2.074, 2.069, 2.064, 2.06, 2.056, 2.052, 2.048, 2.045, 2.042, 1.960};

    public void genMatrix_x(int n, double[][] mainMatrix) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                switch (j) {
                    case 0:
                        if (mainMatrix[i][j] < 0) {
                            matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_min1;
                        } else matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_max1;
                        break;
                    case 1:
                        if (mainMatrix[i][j] < 0) {
                            matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_min2;
                        } else matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_max2;
                        break;
                    case 2:
                        if (mainMatrix[i][j] < 0) {
                            matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_min3;
                        } else matrix_x[i][j] = Math.abs(mainMatrix[i][j]) * x_max3;
                        break;
                    default:
                        break;
                }
            }

        }
        for (int i = 0; i < 15; i++) {
            matrix_x[i][3] = matrix_x[i][0] * matrix_x[i][1];
            matrix_x[i][4] = matrix_x[i][0] * matrix_x[i][2];
            matrix_x[i][5] = matrix_x[i][1] * matrix_x[i][2];
            matrix_x[i][6] = matrix_x[i][0] * matrix_x[i][2] * matrix_x[i][1];
            matrix_x[i][7] = Math.pow(matrix_x[i][0], 2);
            matrix_x[i][8] = Math.pow(matrix_x[i][1], 2);
            matrix_x[i][9] = Math.pow(matrix_x[i][2], 2);
        }
    }

    public double[] countAver_y(double[][] matrix_y) {
        double[] averY = new double[matrix_y.length];
        for (int i = 0; i < matrix_y.length; i++) {
            double currentSum = 0;
            for (int j = 0; j < matrix_y[0].length; j++) {
                currentSum += matrix_y[i][j];
            }
            averY[i] = currentSum / matrix_y[0].length;
        }
        return averY;
    }

    public double[] countMathExp_x(double[][] matrix_x) {
        double[] mathExp_x = new double[matrix_x.length];
        for (int i = 0; i < matrix_x.length; i++) {
            double currentSum = 0;
            for (int j = 0; j < matrix_x[0].length; j++) {
                currentSum += matrix_x[i][j];
            }
            mathExp_x[i] = currentSum / matrix_x[0].length;
        }
        return mathExp_x;
    }

    public double countMathExp_y(double[] averY) {
        double currentSum = 0;
        for (int i = 0; i < averY.length; i++) {
            currentSum += averY[i];
        }
        return currentSum / averY.length;
    }

    public double[] countA_1(double[][] matrix_x, double[] aver_y) {
        double[] A_1 = new double[matrix_x[0].length];
        for (int i = 0; i < matrix_x[0].length; i++) {
            double currentSum = 0;
            for (int j = 0; j < matrix_x.length; j++) {
                currentSum += matrix_x[j][i] * aver_y[j];
            }
            A_1[i] = currentSum / aver_y.length;
        }
        return A_1;
    }


    public void findCoef() {
        double[][] akp = new double[15][10];

        for (int i = 0; i < 15; i++) {
            for (int k = 0; k < 10; k++) {
                for (int p = 0; p < 10; p++) {
                    akp[k][p] += (matrix_x[i][k] * matrix_x[i][p]) / n;
                }
            }
        }
    }

    public double[] calcCoef(double[] mx, double[] a, double mathExp_y) {

        double[][] akp = new double[15][10];

        for (int i = 0; i < 15; i++) {
            for (int k = 0; k < 10; k++) {
                for (int p = 0; p < 10; p++) {
                    akp[k][p] += (matrix_x[i][k] * matrix_x[i][p]) / n;
                }
            }
        }
        double[] arrForAi = {mathExp_y, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9]};

        double[] coef = new double[11];
        double[][] lhsArr = {{1, mx[1], mx[2], mx[3], mx[3], mx[4], mx[5], mx[6], mx[7], 0, 0},
                {mx[1], akp[0][0], akp[1][0], akp[2][0], akp[3][0], akp[4][0], akp[5][0], akp[6][0], akp[7][0], akp[8][0], akp[9][0]},
                {mx[2], akp[0][1], akp[1][1], akp[2][1], akp[3][1], akp[4][1], akp[5][1], akp[6][1], akp[7][1], akp[8][1], akp[9][1]},
                {mx[3], akp[0][2], akp[1][2], akp[2][2], akp[3][2], akp[4][2], akp[5][2], akp[6][2], akp[7][2], akp[8][2], akp[9][2]},
                {mx[4], akp[0][3], akp[1][3], akp[2][3], akp[3][3], akp[4][3], akp[5][3], akp[6][3], akp[7][3], akp[8][3], akp[9][3]},
                {mx[5], akp[0][4], akp[1][4], akp[2][4], akp[3][4], akp[4][4], akp[5][4], akp[6][4], akp[7][4], akp[8][4], akp[9][4]},
                {mx[6], akp[0][5], akp[1][5], akp[2][5], akp[3][5], akp[4][5], akp[5][5], akp[6][5], akp[7][5], akp[8][5], akp[9][5]},
                {mx[7], akp[0][6], akp[1][6], akp[2][6], akp[3][6], akp[4][6], akp[5][6], akp[6][6], akp[7][6], akp[8][6], akp[9][6]},
                {mx[8], akp[0][7], akp[1][7], akp[2][7], akp[3][7], akp[4][7], akp[5][7], akp[6][7], akp[7][7], akp[8][7], akp[9][7]},
                {mx[9], akp[0][8], akp[1][8], akp[2][8], akp[3][8], akp[4][8], akp[5][8], akp[6][8], akp[7][8], akp[8][8], akp[9][8]},
                {mx[10], akp[0][9], akp[1][9], akp[2][9], akp[3][9], akp[4][9], akp[5][9], akp[6][9], akp[7][9], akp[8][9], akp[9][9]}};

        Matrix lhs = new Matrix(lhsArr);
        Matrix rhs = new Matrix(arrForAi, 11);
        //Calculate Solved Matrix
        Matrix ans = lhs.solve(rhs);
        //Printing Answers
        for (int i = 0; i < coef.length - 1; i++) {
            coef[i] = (float) ans.get(i, 0);
        }
        return coef;
    }

    public double[] calcDispersion(double[][] matrix, double[] averMatrix) {
        double currValue;
        double[] Dispersion = new double[averMatrix.length];
        for (int i = 0; i < matrix.length; i++) {
            currValue = 0;
            for (int j = 0; j < matrix[0].length; j++) {
                currValue += Math.pow((matrix[i][j] - averMatrix[i]), 2);
            }
            Dispersion[i] = currValue / matrix[0].length;
        }
        return Dispersion;
    }

    public double criterionCochrane(double[] dispersion) {

        return matrixFunc.max(dispersion) / matrixFunc.sum(dispersion);
    }

    public double[] studentCriterion(int n, double[] dispersion, double[] averY, double[][] coded_x) {
        double[] beta = new double[10];
        double Sb_1 = matrixFunc.sum(dispersion) / n;
        Sb_2 = Sb_1 / (n * m);
        double Sb_3 = Math.sqrt(Sb_2);
        double current;
        for (int i = 0; i < beta.length; i++) {
            current = 0;
            for (int j = 0; j < n; j++) {
                current += coded_x[j][i] * averY[j];
            }
            beta[i] = current / beta.length;
        }
        double[] T = new double[beta.length];
        for (int i = 0; i < beta.length; i++) {
            T[i] = Math.abs(beta[i]) / Sb_3;
        }
        return T;
    }


    double[][] FisheraTable = {
            {164.4, 199.5, 215.7, 224.6, 230.2, 234, 244.9, 249.0, 254.3},
            {18.5, 19.2, 19.2, 19.3, 19.3, 19.3, 19.4, 19.4, 19.5},
            {10.1, 9.6, 9.3, 9.1, 9, 8.9, 8.7, 8.6, 8.5},
            {7.7, 6.9, 6.6, 6.4, 6.3, 6.2, 5.9, 5.8, 5.6},
            {6.6, 5.8, 5.4, 5.2, 5.1, 5, 4.7, 4.5, 4.4},
            {6, 5.1, 4.8, 4.5, 4.4, 4.3, 4, 3.8, 3.7},
            {5.5, 4.7, 4.4, 4.1, 4, 3.9, 3.6, 3.4, 3.2},
            {5.3, 4.5, 4.1, 3.8, 3.7, 3.6, 3.3, 3.1, 2.9},
            {5.1, 4.3, 3.9, 3.6, 3.5, 3.4, 3.1, 2.9, 2.7},
            {5, 4.1, 3.7, 3.5, 3.3, 3.2, 2.9, 2.7, 2.5},
            {4.8, 4, 3.6, 3.4, 3.2, 3.1, 2.8, 2.6, 2.4},
            {4.8, 3.9, 3.5, 3.3, 3.1, 3, 2.7, 2.5, 2.3},
            {4.7, 3.8, 3.4, 3.2, 3, 2.9, 2.6, 2.4, 2.2},
            {4.6, 3.7, 3.3, 3.1, 3, 2.9, 2.5, 2.3, 2.1},
            {4.5, 3.7, 3.3, 3.1, 2.9, 2.8, 2.5, 2.3, 2.1},
            {4.5, 3.6, 3.2, 3, 2.9, 2.7, 2.4, 2.2, 2},
            {4.5, 3.6, 3.2, 3, 2.8, 2.7, 2.4, 2.2, 2},
            {4.4, 3.6, 3.2, 2.9, 2.8, 2.7, 2.3, 2.1, 1.9},
            {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1, 1.9},
            {4.4, 3.5, 3.1, 2.9, 2.7, 2.6, 2.3, 2.1, 1.9},
            {4.3512, 3.4928, 3.0984, 2.8661, 2.7109, 2.5990, 2.5140, 2.4471, 1.8},
            {4.3248, 3.4668, 3.0725, 2.8401, 2.6848, 2.5727, 2.4876, 2.4205, 1.8},
            {4.3009, 3.4434, 3.0491, 2.8167, 2.6613, 2.5491, 2.4638, 2.3965, 1.7},
            {4.2793, 3.4221, 3.0280, 2.7955, 2.6400, 2.7763, 2.6207, 2.5082, 1.7},
            {4.2597, 3.4028, 3.0088, 2.7763, 2.6207, 2.5082, 2.4226, 2.3551, 1.7},
            {4.2417, 3.3852, 2.9912, 2.7587, 2.6030, 2.4904, 2.4047, 2.3371, 1.7},
            {4.2252, 3.3690, 2.9752, 2.7426, 2.5868, 2.4741, 2.3883, 2.3205, 1.7},
            {4.2100, 3.3541, 2.9604, 2.7278, 2.5719, 2.4591, 2.3732, 2.3053, 1.7},
            {4.1960, 3.3404, 2.9467, 2.7141, 2.5581, 2.4453, 2.3593, 2.2913, 1.6},
            {4.1830, 3.3277, 2.9340, 2.7014, 2.5454, 2.4324, 2.3463, 2.2783, 1.6},
            {4.1709, 3.3158, 2.9223, 2.6896, 2.5336, 2.4205, 2.3343, 2.2662, 1.5},
            {3.8, 3.0, 2.6, 2.4, 2.2, 2.1, 1.8, 1.5, 1},
    };

    public double fisherCriterion(double[] T, double[][] array_x, double[] coef, double[] averY) {
        double[] y_arr_1 = new double[15];
        double[][] currentArr = array_x;
        double currentSum;
        int d = 0;
//      Ця перевірка використанна для того щоб порахувати кількість значущих коефіцієнтів
        for (int i = 0; i < T.length; i++) {
            if (T[i] > StudentaTable[m - 1]) {
                d++;
            }
        }
        System.out.println("d:" + d);
        for (int i = 1; i < currentArr[0].length; i++) {
            for (int j = 0; j < currentArr[0].length; j++) {
                if (T[i] > StudentaTable[m - 1]) {                         //Для перевірки ми використовуємо лише значущі коефіцієнти регресії
                    currentArr[i][j] *= coef[i];
                }
            }
        }
        for (int i = 0; i < currentArr[0].length; i++) {
            currentSum = 0;
            for (int j = 0; j < currentArr[0].length; j++) {
                currentSum += currentArr[i][j];
            }
            if (T[0] > StudentaTable[m - 1]) {                            //Для перевірки ми використовуємо лише значущі коефіцієнти регресії
                currentSum += coef[0];
            }
            y_arr_1[i] = currentSum;

        }
        double Sad = 0;
        for (int i = 0; i < averY.length - 1; i++) {
            Sad += Math.pow(y_arr_1[i] - averY[i], 2);
        }
        Sad = Sad * m / (n - d);
//        System.out.println(Sb_2);
        return Sad / (Sb_2 * constanta);
    }
}