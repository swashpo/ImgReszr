/* *****************************************************************************
 *  Name:    Richard Ma
 *  NetID:   rm44
 *  Precept: P05
 *
 *
 *  Description:  Given a picture, finds horizontal and vertical seams with
 *  lowest energy and also removes specified seams from picture.
 *
 **************************************************************************** */

import edu.princeton.cs.algs4.IndexMinPQ;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Stopwatch;

public class SeamCarver {
    private Picture pic; // photo

    private int width; // of pict
    private int height; // of pict

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException("No pic");

        pic = new Picture(picture);
        width = pic.width();
        height = pic.height();
    }

    // current picture
    public Picture picture() {
        return new Picture(pic);
    }

    // width of current picture
    public int width() {
        return width;
    }

    // height of current picture
    public int height() {
        return height;
    }

    // computes the delta value given two pixel colors
    private double delta(int x, int y) {
        // compute rgb of x and y
        int r = (x >> 16) & 0xFF;
        int g = (x >> 8) & 0xFF;
        int b = (x >> 0) & 0xFF;

        int r2 = (y >> 16) & 0xFF;
        int g2 = (y >> 8) & 0xFF;
        int b2 = (y >> 0) & 0xFF;

        double delta = Math.pow(r - r2, 2) + Math.pow(g - g2, 2) +
                Math.pow(b - b2, 2);

        return delta;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= width)
            throw new IllegalArgumentException("Outside range");
        if (y < 0 || y >= height)
            throw new IllegalArgumentException("Outside range");

        // left pixel
        int lt = (x - 1 + width) % width;
        // right pixel
        int rt = (x + 1) % width;

        lt = pic.getRGB(lt, y);
        rt = pic.getRGB(rt, y);

        double deltaX = delta(lt, rt);

        // top pixel
        int tp = (y - 1 + height) % height;
        // bottom pixel
        int bt = (y + 1) % height;

        tp = pic.getRGB(x, tp);
        bt = pic.getRGB(x, bt);

        double deltaY = delta(tp, bt);

        double energy = Math.sqrt(deltaX + deltaY);
        return energy;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        transpose();
        int[] seam = findVerticalSeam();
        transpose();

        return seam;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        // stores pixels to be examined
        IndexMinPQ<Double> least = new IndexMinPQ<Double>(width * height);
        double[][] energyTable = new double[height][width];

        // calculating energies in picture
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                double energy = energy(i, j);
                energyTable[j][i] = energy;
            }
        }

        // array to keep track of what leads to what
        int[][] pathTo = new int[height][width];
        for (int i = 0; i < width; i++) {
            pathTo[0][i] = i;
        }

        // array to keep track of energy totals to each cell
        double[][] energyTo = new double[height][width];

        // set first row energy to
        for (int i = 0; i < width; i++) {
            energyTo[0][i] = energyTable[0][i];
        }

        // set remaining to infinity
        for (int i = 1; i < height; i++) {
            for (int j = 0; j < width; j++) {
                energyTo[i][j] = Double.POSITIVE_INFINITY;
            }
        }

        // we start with top virtual root by inserting all top row
        for (int i = 0; i < width; i++) {
            least.insert(i, energyTable[0][i]);
        }

        while (!least.isEmpty()) {
            int lowest = least.delMin();
            int row = lowest / width;
            int col = lowest % width;

            // update total energies for each neighboring pixel

            // directly below
            if ((row + 1) < height) {
                relax(row, col, 0, energyTo, pathTo, lowest, energyTable, least);
            }

            // below and to the left
            if ((row + 1) < height && col > 0) {
                relax(row, col, -1, energyTo, pathTo, lowest, energyTable, least);
            }

            // below and to the right
            if ((row + 1) < height && col + 1 < width) {
                relax(row, col, 1, energyTo, pathTo, lowest, energyTable, least);
            }
        }

        // compute path

        // find which pixel in final row as lowest energyTo
        double low = Double.POSITIVE_INFINITY;
        int lowRow = 0;
        for (int i = 0; i < width; i++) {
            if (energyTo[height - 1][i] < low) {
                low = energyTo[height - 1][i];
                lowRow = i;
            }
        }

        int[] path = new int[height]; // path to be carved
        path[height - 1] = lowRow; // last in array

        // trace back from the bottom w lowest energyTo
        // store all points in array of size 'height' from last to first index
        // return array

        for (int i = height - 2; i > -1; i--) {
            int marker = pathTo[i + 1][lowRow];
            lowRow = marker % width;
            path[i] = lowRow;
        }
        return path;
    }

    // relax or go to neighbor
    private void relax(int row, int col, int path, double[][] energyTo,
                       int[][] pathTo, int lowest, double[][] energyTable,
                       IndexMinPQ<Double> least) {

        double temp = energyTo[row][col] + energyTable[row + 1][col + path];

        // if updated energy is less, update table
        if (temp < energyTo[row + 1][col + path]) {
            pathTo[row + 1][col + path] = lowest;
            energyTo[row + 1][col + path] = temp;

            int marker = (row + 1) * width + col + path;

            if (least.contains(marker))
                least.decreaseKey(marker, temp);
            else
                least.insert(marker, temp);
        }
    }

    // transposes photo
    private void transpose() {
        int newHt = width;
        int newWt = height;

        // new pic with flipped dim
        Picture transpose = new Picture(newWt, newHt);

        for (int i = 0; i < newHt; i++) {
            for (int j = 0; j < newWt; j++) {
                transpose.setRGB(j, i, pic.getRGB(i, j));
            }
        }

        pic = transpose;
        int temp = width;
        width = height;
        height = temp;
    }


    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("No seam");
        if (height == 1) throw new IllegalArgumentException("Can't!");
        if (seam.length != width)
            throw new IllegalArgumentException("Not a seam");

        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= height)
                throw new IllegalArgumentException("Not valid");
            if (i == seam.length - 1) continue;
            else if (Math.abs(seam[i] - seam[i + 1]) > 1)
                throw new IllegalArgumentException("Not valid");
        }

        transpose();
        removeVerticalSeam(seam);
        transpose();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("No seam");
        if (width == 1) throw new IllegalArgumentException("Can't!");
        if (seam.length != height)
            throw new IllegalArgumentException("Not a seam");

        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= width)
                throw new IllegalArgumentException("Not valid");
            if (i == seam.length - 1) continue;
            else if (Math.abs(seam[i] - seam[i + 1]) > 1)
                throw new IllegalArgumentException("Not valid");
        }

        width--;
        Picture update = new Picture(width, height);

        // create new pic object with one less width dim
        // add in all pixels, skipping over ones to be removed

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (j < seam[i]) {
                    update.setRGB(j, i, pic.getRGB(j, i));
                }
                else
                    update.setRGB(j, i, pic.getRGB(j + 1, i));
            }
        }
        pic = update;
    }

    // main
    public static void main(String[] args) {

        if (args.length != 3) {
            StdOut.println(
                    "Usage:\njava ResizeDemo [image filename] [num columns to remove] [num rows to remove]");
            return;
        }

        Picture picture = new Picture(args[0]);
        int removeColumns = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]);

        StdOut.printf("%d-by-%d image\n", picture.width(), picture.height());
        SeamCarver sc = new SeamCarver(picture);

        Stopwatch sw = new Stopwatch();

        for (int i = 0; i < removeRows; i++) {
            int[] horizontalSeam = sc.findHorizontalSeam();
            sc.removeHorizontalSeam(horizontalSeam);
        }

        for (int i = 0; i < removeColumns; i++) {
            int[] verticalSeam = sc.findVerticalSeam();
            sc.removeVerticalSeam(verticalSeam);
        }

        StdOut.printf("new image size is %d columns by %d rows\n", sc.width(), sc.height());

        StdOut.println("Resizing time: " + sw.elapsedTime() + " seconds.");
        picture.show();
        sc.picture().show();

        StdOut.println(sc.energy(5, 5));
    }
}

