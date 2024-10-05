import javax.swing.*;
import java.awt.*;
import java.awt.image.AreaAveragingScaleFilter;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;
import java.util.zip.CheckedOutputStream;

/**
 * Program to calculate and draw the shortest distance to travel
 * to all the coordinates and return to starting point.
 * Coordinates are taken from a file.
 * Program uses 2 algorithms: Brute-Force and Ant Colony Optimization.
 * To use Brute-Force, set variable "chosenMethod = 1"
 * To use Ant Colony Optimization, set variable "chosenMethod = 2"
 * @author Yeldos Urashev
 * @since 12.05.2024
 */

public class YeldosUrashev {

    public static final double alpha = 0.6;
    public static final int N = 100;
    public static final int M = 50;
    public static final double degrade = 0.9;
    public static final double beta = 10;
    public static double initialPheromone = 0.1;
    public static final double Q = 0.0001;

    public static void main(String[] args) throws FileNotFoundException {

        int chosenMethod = 1; // Choose 1 to brute-force and 2 for ACO
        int choiceOfDrawing = 1; // Choose 1 to draw the shortest path and 2 to draw pheromone graph for ACO
        String fileName = "input01.txt"; // Change the input file by changing this variable

        int houseNumber = 0;
        double[] minDistance = {Double.MAX_VALUE};
        int width = 900, height = 900;
        double penRadius = 0.005, circleRadius = 0.02;
        ArrayList<Integer> finalRoute = new ArrayList<>();
        double[][] pheromones = new double[0][];

        ArrayList<Double[]> houses = new ArrayList<>();

        /*
        Getting the coordinates from a file
        Checking if it was found
         */
        File coordinates = new File(fileName);
        if (!coordinates.exists()) {
            System.out.printf("%s cannot be found", coordinates);
            System.exit(1);
        }
        /*
        Scanning the coordinates file
        Saving all x, y values to ArrayList
         */
        Scanner scanCoordinates = new Scanner(coordinates);

        while (scanCoordinates.hasNextLine()) {
        String[] line = scanCoordinates.nextLine().split(",");
        double x = Double.parseDouble(line[0]);
        double y = Double.parseDouble(line[1]);
        Double[] temp = new Double[2];
        temp[0] = x;
        temp[1] = y;
        houses.add(temp);
        houseNumber += 1;
        }
        scanCoordinates.close();

        long start = System.currentTimeMillis(); // starting to count the time of compilation
        houseNumber--;
        int[] routes = new int[houseNumber];
        int[] index = new int[houseNumber];
        for (int i = 0; i < houseNumber; i++) {
            index[i] = i;
        }


        if(chosenMethod==1) {
            /*
            Brute-force method to compute all the possible routes recursively.
            */
            permute(routes, index, 0, houseNumber-1, houses, minDistance);
            finalRoute.add(1);
            for(int i:routes) {
                finalRoute.add(i+2);
            }
            finalRoute.add(1);
        }
        else {
            /*
            Implementation of ACO algorithm
             */
            houseNumber++;
            double tempDistance = 0;
            double delta = 0;
            double[][] edgeValues = new double[houseNumber][houseNumber];
            pheromones = new double[houseNumber][houseNumber];
            int[] route = new int[houseNumber];
            for(int i=0; i<houseNumber; i++) {
                route[i] = i;
            }
            Random rand = new Random();
            setInitialPheromone(pheromones, initialPheromone); // setting initial pheromones
            setInitialEdges(edgeValues, houses, pheromones); // setting initial edge values

            for (int i = 0; i < N; i++) { // Loop for iterations
                for (int j = 0; j < M; j++) { // Loop for ants
                    updateEdges(edgeValues, houses, pheromones, route); // after every cycle, we need to update the edge values

                    int startNode = rand.nextInt(houseNumber); // putting an ant to a random coordinate
                    tempDistance = cycle(startNode, edgeValues, route, houses); // getting the total distance for current cycle
                    if(tempDistance < minDistance[0]) {
                        minDistance[0] = tempDistance; // updating the shortest distance found
                        routes = route.clone();
                    }

                    delta = Q / tempDistance;
                    updatePheromones(pheromones, delta, route); // after every cycle, we need to update pheromone values
                    //show(route);
                }
                degradePheromones(pheromones); // after every iteration, we need to simulate evaporation of past pheromones
            }

            /*
            Reordering the shortest route.
             */
            int startIndex = 0;
            for (int i=0; i<routes.length; i++) {
                if(routes[i]==0) startIndex = i;
            }
            for(int i=startIndex; i<routes.length; i++) {
                finalRoute.add(routes[i]+1);
            }
            for(int i=0; i<startIndex; i++) {
                finalRoute.add(routes[i]+1);
            }
            finalRoute.add(1);
        }

        long end = System.currentTimeMillis(); // ending the calculation of compilation time


        StdDraw.setCanvasSize(width, height);
        StdDraw.setXscale(0, 1.0);
        StdDraw.setYscale(0, 1.0);
        StdDraw.enableDoubleBuffering();
        if(chosenMethod==1 || (chosenMethod==2 && choiceOfDrawing==1)) {
            /*
            Drawing the shortest route
             */
            StdDraw.setPenRadius(penRadius);
            Color color = Color.LIGHT_GRAY;
            int prev = 0;
            StdDraw.setPenColor(Color.BLACK);
            for(int i:finalRoute) {
                StdDraw.line(houses.get(i-1)[0], houses.get(i-1)[1], houses.get(prev)[0], houses.get(prev)[1]);
                prev = i-1;
            }
            for(int i:finalRoute) {
                if(i==1) color = StdDraw.PRINCETON_ORANGE;
                else color = StdDraw.LIGHT_GRAY;
                StdDraw.setPenColor(color);
                StdDraw.filledCircle(houses.get(i-1)[0], houses.get(i-1)[1], circleRadius);
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(houses.get(i-1)[0], houses.get(i-1)[1], String.valueOf(i));
            }
            StdDraw.show();
            if(chosenMethod==1) {
                System.out.println("Method: Brute-Force Method");
            }
            else if(chosenMethod==2) {
                System.out.println("Method: Ant Colony Optimization Method");
            }
            System.out.println("Shortest Distance: " + minDistance[0]);
            System.out.print("Shortest path: ");
            System.out.print(finalRoute);
            System.out.println();
            System.out.println("Time it takes to find the shortest path: " + (end-start)/1000.0 + " seconds.");
        }
        else{
            /*
            Drawing the pheromone intensity graph
             */
            Color color = Color.LIGHT_GRAY;
            StdDraw.setPenColor(Color.BLACK);
            for(int i=0; i < pheromones.length; i++) {
                for(int j=0; j < pheromones[i].length; j++) {
                    StdDraw.setPenRadius(pheromones[i][j]/1.2);
                    StdDraw.line(houses.get(i)[0], houses.get(i)[1], houses.get(j)[0], houses.get(j)[1]);
                }
            }
            StdDraw.setPenRadius(penRadius);
            for(int i:finalRoute) {
                if(i==1) color = StdDraw.PRINCETON_ORANGE;
                else color = StdDraw.LIGHT_GRAY;
                StdDraw.setPenColor(color);
                StdDraw.filledCircle(houses.get(i-1)[0], houses.get(i-1)[1], circleRadius);
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(houses.get(i-1)[0], houses.get(i-1)[1], String.valueOf(i));
            }
            StdDraw.show();
            System.out.println("Method: Ant Colony Optimization Method");
            System.out.println("Shortest Distance: " + minDistance[0]);
            System.out.print("Shortest path: ");
            System.out.print(finalRoute);
            System.out.println();
            System.out.println("Time it takes to find the shortest path: " + (end-start)/1000.0 + " seconds.");

        }

    }

    /**
     * Sets initial pheromone values
     * @param pheromones array to keep pheromone values
     * @param initialPheromone initial pheromone hyperparameter
     */
    private static void setInitialPheromone(double[][] pheromones, double initialPheromone) {
        for(int i = 0; i < pheromones.length; i++) {
            for(int j=i+1; j < pheromones.length; j++) {
                pheromones[i][j] = initialPheromone;
                pheromones[j][i] = initialPheromone;
            }
        }
    }

    /**
     * Sets initial edge values
     * @param edgeValues array to keep edge values
     * @param houses array with all the coordinates of houses
     * @param pheromones array to keep pheromone values
     */
    private static void setInitialEdges(double[][] edgeValues, ArrayList<Double[]> houses, double[][] pheromones) {
        for(int i=0; i < edgeValues.length; i++) {
            for(int j=i+1; j < edgeValues.length; j++) {
                double distance = calculateDistance(houses.get(i), houses.get(j));
                edgeValues[i][j] =  Math.pow(pheromones[i][j], alpha) / Math.pow(distance, beta);
                edgeValues[j][i] =  Math.pow(pheromones[j][i], alpha) / Math.pow(distance, beta);
            }
        }
    }

    /**
     * Decreases pheromone values by some factor
     * @param pheromones array to keep pheromone values
     */
    private static void degradePheromones(double[][] pheromones) {
        for (int i = 0; i < pheromones.length; i++) {
            for (int j = i + 1; j < pheromones[i].length; j++) {
                pheromones[i][j] = pheromones[i][j] * degrade;
                pheromones[j][i] = pheromones[j][i] * degrade;
            }
        }
    }

    /**
     * Simulates a cycle of one ant
     * @param startNode the index where ant starts
     * @param edgeValues array to keep edge values
     * @param route array to keep the route of an ant
     * @param houses array with all the coordinates of houses
     * @return the total distance traveled by an ant in this cycle
     */
    private static double cycle(int startNode, double[][] edgeValues, int[] route, ArrayList<Double[]> houses) {
        int firstNode = startNode;
        double[] probabilities = new double[edgeValues.length];
        boolean[] visited = new boolean[edgeValues.length];
        int j = 0;
        route[j++] = startNode;
        visited[startNode] = true;
        int nextNode = 0;
        double distance = 0;

        for(int i=0; i < edgeValues.length-1; i++) {
            getProbabilities(probabilities, visited, edgeValues, startNode);
            nextNode = getNextNode(probabilities);
            visited[nextNode] = true;
            route[j++] = nextNode;
            distance += calculateDistance(houses.get(startNode), houses.get(nextNode));
            startNode = nextNode;
        }
        distance += calculateDistance(houses.get(firstNode), houses.get(nextNode));
        return distance;
    }

    /**
     * Calculates probabilistically the next node to go
     * @param probabilities array with the probabilities of all edges
     * @return the index of the next node or -1 if not found
     */
    private static int getNextNode(double[] probabilities) {
        double probability = Math.random();
        double cumulativeProbability = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulativeProbability += probabilities[i];
            if(probability <= cumulativeProbability) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Calculates the probabilities of all the nodes, by their edge distances
     * @param probabilities array to keep the probabilities
     * @param visited array keeping visited nodes
     * @param edgeValues array keeping edge values
     * @param startNode the node where an ant is located
     */
    private static void getProbabilities(double[] probabilities, boolean[] visited, double[][] edgeValues, int startNode) {
        double edgeSum = getEdgeSum(visited, edgeValues, startNode);
        for(int i = 0; i < probabilities.length; i++) {
            if(visited[i]) probabilities[i] = 0;
            else probabilities[i] = edgeValues[startNode][i] / edgeSum;
        }

    }


    /**
     * Calculates total sum of edges that connect starting node to unvisited nodes
     * @param visited array keeping visited nodes
     * @param edgeValues array keeping edge values
     * @param startNode the node where an ant is located
     * @return sum of edge values
     */
    private static double getEdgeSum(boolean[] visited, double[][] edgeValues, int startNode) {
        double sum = 0;
        for(int i=0; i < visited.length; i++) {
            if(!visited[i]) sum += edgeValues[startNode][i];
        }
        return sum;
    }

    /**
     * Updates pheromone values using relevant delta value
     * @param pheromones array to keep pheromone values
     * @param delta constant that is calculated by formula Q / totalCycleDistance
     * @param route array keeping the relevant route of an ant
     */
    private static void updatePheromones(double[][] pheromones, double delta, int[] route) {
        int initial = route[0];
        for(int i=1; i<route.length; i++) {
            pheromones[initial][route[i]] += delta;
            pheromones[route[i]][initial] += delta;
            initial = route[i];
        }

    }

    /**
     * Updates edge values using relevant pheromone values
     * @param edgeValues array to keep the edge values
     * @param houses array keeping the coordinates of houses
     * @param pheromones array keeping the relevant pheromone values
     * @param route array keeping the relevant route of an ant
     */
    private static void updateEdges(double[][] edgeValues, ArrayList<Double[]> houses, double[][] pheromones, int[] route) {
        int initial = route[0];
        for(int i=1; i<route.length; i++) {
            double distance = calculateDistance(houses.get(initial), houses.get(route[i]));
            edgeValues[initial][route[i]] =  Math.pow(pheromones[initial][route[i]], alpha) / Math.pow(distance, beta);
            edgeValues[route[i]][initial] =  Math.pow(pheromones[route[i]][initial], alpha) / Math.pow(distance, beta);
        }
    }

    /**
     * Goes through all the permutations and keeps the shortest distance and route
     * @param routes array to keep the shortest route
     * @param index default array to permute
     * @param start starting index of an array to permute
     * @param end last index of an array to permute
     * @param houses array keeping the coordinates of houses
     * @param minDistance variable to keep the shortest distance found
     */
    private static void permute(int[] routes, int[] index, int start, int end, ArrayList<Double[]> houses, double[] minDistance) {
        if (start == end) {
            double distance = calculateTotalDistance(index, houses);
            if (minDistance[0] > distance) {
                minDistance[0] = distance;
                for(int i=0; i<index.length; i++) {
                    routes[i] = index[i];
                }
            }
            return;
        }

        for (int i = start; i <= end; i++) {
            int temp = index[start];
            index[start] = index[i];
            index[i] = temp;
            permute(routes, index, start + 1, end, houses, minDistance);
            temp = index[start];
            index[start] = index[i];
            index[i] = temp;
        }
    }

    /**
     * Calculates total distance with an array of coordinates
     * @param route array keeping the indexes in route order
     * @param houses array keeping all the coordinates of houses
     * @return the total distance
     */
    public static double calculateTotalDistance(int[] route, ArrayList<Double[]> houses) {
        double distance = 0, x = 0, y = 0;
        int prev = -1;
        for (int id : route) {
            x = Math.pow(houses.get(prev+1)[0] - houses.get(id+1)[0], 2);
            y = Math.pow(houses.get(prev+1)[1] - houses.get(id+1)[1], 2);
            distance += Math.pow(x+y, 0.5);
            prev = id;
        }
        x = Math.pow(houses.get(0)[0] - houses.get(prev+1)[0], 2);
        y = Math.pow(houses.get(0)[1] - houses.get(prev+1)[1], 2);
        distance += Math.pow(x+y, 0.5);

        return distance;
    }

    /**
     * Calculates the distance between two nodes
     * @param first x,y coordinate of first node
     * @param second x,y coordinate of second node
     * @return distance between two nodes
     */
    public static double calculateDistance(Double[] first, Double[] second) {
        double distance = 0;
        double x = first[0] - second[0];
        double y = first[1] - second[1];
        distance = Math.pow(x*x + y*y, 0.5);
        return distance;
    }
}