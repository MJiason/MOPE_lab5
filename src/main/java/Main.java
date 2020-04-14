public class Main {
    public static void main(String[] args) {
        int n = 15;
        int m = 3;
        Experiment experiment = new Experiment(-2,-10, -3,4,8,6);
        MatrixFunc matrixFunc = new MatrixFunc();
        experiment.genMatrix_x(15, experiment.mainMatrix);
        double[][] matrixY =  experiment.genMatrix();
        double[] averY = experiment.countAver_y(matrixY);
        double[] averX = experiment.countAver_y(experiment.matrix_x);
        double mathExp_y = experiment.countMathExp_y(averY);
        double[] A_1 = experiment.countA_1(experiment.matrix_x, averY);
        double[] coef = experiment.calcCoef(averX, A_1, mathExp_y);
        double[] dispersion = experiment.calcDispersion(matrixY, averY);
        double criterionCochrane = experiment.criterionCochrane(dispersion);
        double[] criterionStudent = experiment.studentCriterion(n, dispersion, averY, experiment.mainMatrix);
        double fisherCriterion = experiment.fisherCriterion(criterionStudent, experiment.matrix_x,coef, averY);
        matrixFunc.printMatrix("Matrix_x", experiment.matrix_x);
        matrixFunc.printMatrix("Matrix_y", matrixY);
        matrixFunc.printMatrix("averY", averY);
        matrixFunc.printMatrix("averX", averX);
        System.out.format("MathExpectationY %2.2f\n", mathExp_y);
        matrixFunc.printMatrix("A_1", A_1);
        matrixFunc.printMatrix("Coeficients", coef);
        matrixFunc.printMatrix("Dispersion", dispersion);
        System.out.format("criterionCochrane:\nf1: " + (m-1) + "\nf2: " + (n-1) + "\nGp: %2.2f, отже дисперсія однорідна", criterionCochrane);
        System.out.format("\nStudentCriterion:\nf3=f1*f2 =(m-1)*N: " + 30);
        matrixFunc.printMatrix("\nt", criterionStudent);
        System.out.println("FisherCriterion:\nf4 = N – d=6\nf3 = f1*f2=(m-1)*N=2*6=12: " + fisherCriterion);

    }

}
