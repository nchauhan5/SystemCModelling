/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 5 solution */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "systemc.h"

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */
#define SET_STACK_SIZE set_stack_size(128*1024*1024);
/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

typedef struct Image_s {
	unsigned char img[SIZE];

	Image_s(void) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = 0;
		}
	}

	Image_s& operator=(const Image_s& copy) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = copy.img[i];
		}
		return *this;
	}

	operator unsigned char*() {
		return img;
	}

	unsigned char& operator[](const int index) {
		return img[index];
	}
} IMAGE;

typedef struct SImage_s {
	short int img[SIZE];

	SImage_s(void) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = 0;
		}
	}

	SImage_s& operator=(const SImage_s& copy) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = copy.img[i];
		}
		return *this;
	}

	operator short int*() {
		return img;
	}

	short int& operator[](const int index) {
		return img[index];
	}
} SIMAGE;

typedef struct KERNEL_s {
	float kernel[WINSIZE];

	KERNEL_s(void) {
		for (int i = 0; i < WINSIZE; i++) {
			kernel[i] = 0;
		}
	}

	KERNEL_s& operator=(const KERNEL_s& copy) {
		for (int i = 0; i < WINSIZE; i++) {
			kernel[i] = copy.kernel[i];
		}
		return *this;
	}

	operator float*() {
		return kernel;
	}

	float& operator[](const int index) {
		return kernel[index];
	}
} KERNEL;

typedef struct TEMPIM_s {
	float tempim[SIZE];

	TEMPIM_s(void) {
		for (int i = 0; i < SIZE; i++) {
			tempim[i] = 0;
		}
	}

	TEMPIM_s& operator=(const TEMPIM_s& copy) {
		for (int i = 0; i < SIZE; i++) {
			tempim[i] = copy.tempim[i];
		}
		return *this;
	}

	operator float*() {
		return tempim;
	}

	float& operator[](const int index) {
		return tempim[index];
	}
} TEMPIM;

SC_MODULE(Stimulus) {
	IMAGE imageout;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo_out<sc_time> frameTimeOut;

	sc_time startTime;
	/******************************************************************************
	 * Function: read_pgm_image
	 * Purpose: This function reads in an image in PGM format. The image can be
	 * read in from either a file or from standard input. The image is only read
	 * from standard input when infilename = NULL. Because the PGM format includes
	 * the number of columns and the number of rows in the image, these are read
	 * from the file. Memory to store the image is allocated OUTSIDE this function.
	 * The found image size is checked against the expected rows and cols.
	 * All comments in the header are discarded in the process of reading the
	 * image. Upon failure, this function returns 0, upon sucess it returns 1.
	 ******************************************************************************/
	int read_pgm_image(const char *infilename, unsigned char *image, int rows,
			int cols) {
		FILE *fp;
		char buf[71];
		int r, c;

		/***************************************************************************
		 * Open the input image file for reading if a filename was given. If no
		 * filename was provided, set fp to read from standard input.
		 ***************************************************************************/
		if (infilename == NULL)
			fp = stdin;
		else {
			if ((fp = fopen(infilename, "r")) == NULL) {
				fprintf(stderr,
						"Error reading the file %s in read_pgm_image().\n",
						infilename);
				return (0);
			}
		}

		/***************************************************************************
		 * Verify that the image is in PGM format, read in the number of columns
		 * and rows in the image and scan past all of the header information.
		 ***************************************************************************/
		fgets(buf, 70, fp);
		if (strncmp(buf, "P5", 2) != 0) {
			fprintf(stderr, "The file %s is not in PGM format in ", infilename);
			fprintf(stderr, "read_pgm_image().\n");
			if (fp != stdin)
				fclose(fp);
			return (0);
		}
		do {
			fgets(buf, 70, fp);
		} while (buf[0] == '#'); /* skip all comment lines */
		sscanf(buf, "%d %d", &c, &r);
		if (c != cols || r != rows) {
			fprintf(stderr, "The file %s is not a %d by %d image in ",
					infilename, cols, rows);
			fprintf(stderr, "read_pgm_image().\n");
			if (fp != stdin)
				fclose(fp);
			return (0);
		}
		do {
			fgets(buf, 70, fp);
		} while (buf[0] == '#'); /* skip all comment lines */

		/***************************************************************************
		 * Read the image from the file.
		 ***************************************************************************/
		if ((unsigned) rows != fread(image, cols, rows, fp)) {
			fprintf(stderr,
					"Error reading the image data in read_pgm_image().\n");
			if (fp != stdin)
				fclose(fp);
			return (0);
		}

		if (fp != stdin)
			fclose(fp);
		return (1);
	}

	void main(void) {
		// cout << "******In Stimulus**********" << endl;
		int i = 0, n = 0;
		char infilename[40];

		for (i = 0; i < IMG_NUM; i++) {
			n = i % AVAIL_IMG;
			sprintf(infilename, IMG_IN, n + 1);

			/****************************************************************************
			 * Read in the image.
			 ****************************************************************************/
			if (VERBOSE)
				printf("Reading the image %s.\n", infilename);
			if (read_pgm_image(infilename, imageout, ROWS, COLS) == 0) {
				fprintf(stderr, "Error reading the input image, %s.\n",
						infilename);
				exit(1);
			}
			ImgOut.write(imageout);
			startTime = sc_time_stamp();
			int start = startTime.to_seconds() * 1000;
			printf("%d ms : Stimulus sent frame %d.\n", start, i + 1);
			frameTimeOut.write(startTime);
		}
	}

	SC_CTOR(Stimulus) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(Monitor) {
	IMAGE imagein;
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_in<sc_time> frameTimeIn;

	sc_time startTime;
	sc_time endTime;
	/******************************************************************************
	 * Function: write_pgm_image
	 * Purpose: This function writes an image in PGM format. The file is either
	 * written to the file specified by outfilename or to standard output if
	 * outfilename = NULL. A comment can be written to the header if coment != NULL.
	 ******************************************************************************/
	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
			int cols, const char *comment, int maxval) {
		FILE *fp;

		/***************************************************************************
		 * Open the output image file for writing if a filename was given. If no
		 * filename was provided, set fp to write to standard output.
		 ***************************************************************************/
		if (outfilename == NULL)
			fp = stdout;
		else {
			if ((fp = fopen(outfilename, "w")) == NULL) {
				fprintf(stderr,
						"Error writing the file %s in write_pgm_image().\n",
						outfilename);
				return (0);
			}
		}

		/***************************************************************************
		 * Write the header information to the PGM file.
		 ***************************************************************************/
		fprintf(fp, "P5\n%d %d\n", cols, rows);
		if (comment != NULL)
			if (strlen(comment) <= 70)
				fprintf(fp, "# %s\n", comment);
		fprintf(fp, "%d\n", maxval);

		/***************************************************************************
		 * Write the image data to the file.
		 ***************************************************************************/
		if ((unsigned) rows != fwrite(image, cols, rows, fp)) {
			fprintf(stderr,
					"Error writing the image data in write_pgm_image().\n");
			if (fp != stdout)
				fclose(fp);
			return (0);
		}

		if (fp != stdout)
			fclose(fp);
		return (1);
	}

	void main(void) {
		// cout << "******In Monitor**********" << endl;
		char outfilename[128]; /* Name of the output "edge" image */
		int i, n;
		sc_time previousFrame = SC_ZERO_TIME;
		sc_time nextFrame = SC_ZERO_TIME;

		for (i = 0; i < IMG_NUM; i++) {
			frameTimeIn.read(startTime);
			// endTime = sc_time_stamp();
			// int end = endTime.to_seconds() * 1000;
			// int diff = (endTime.to_seconds() - startTime.to_seconds()) * 1000;
			ImgIn.read(imagein);

			/****************************************************************************
			 * Write out the edge image to a file.
			 ****************************************************************************/
			n = i % AVAIL_IMG;
			sprintf(outfilename, IMG_OUT, n + 1);
			if (VERBOSE)
				printf("Writing the edge image in the file %s.\n", outfilename);
			if (write_pgm_image(outfilename, imagein, ROWS, COLS, "", 255)
					== 0) {
				fprintf(stderr, "Error writing the edge image, %s.\n",
						outfilename);
				exit(1);
			}
			cout << sc_time_stamp() << " : Monitor received frame " << i + 1
					<< " with " << sc_time_stamp() - startTime << " delay.\n";
			previousFrame = nextFrame;
			nextFrame = sc_time_stamp();
			if (i >= 1) {
			cout << (nextFrame - previousFrame).to_seconds()
			<< " seconds after previous frame " << i<<" , "
			<< 1 / (nextFrame - previousFrame).to_seconds()
			<< " FPS.\n";
		}
	}
	if (VERBOSE)
		printf("Monitor exits simulation.\n");
	sc_stop();	// done testing, quit the simulation
}

SC_CTOR(Monitor) {
	SC_THREAD(main);
	SET_STACK_SIZE
}
};

SC_MODULE(DataIn) {
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main() {
		// cout << "******In DataIn**********" << endl;
		while (1) {
			ImgIn.read(Image);
			ImgOut.write(Image);
		}
	}

	SC_CTOR(DataIn) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(DataOut) {
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main() {
		// cout << "******In DataOut**********" << endl;
		while (1) {
			ImgIn.read(Image);
			ImgOut.write(Image);
		}
	}

	SC_CTOR(DataOut) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(MAKE_GAUSSIAN_KERNEL) {
	IMAGE imgIn;

	KERNEL ker;
	int winsize;

	sc_fifo_in<IMAGE> imagein;

	sc_fifo_out<KERNEL> kernel;
	sc_fifo_out<int> windowsize;

	sc_fifo_out<IMAGE> imageOut;
	/*******************************************************************************
	 * PROCEDURE: make_gaussian_kernel
	 * PURPOSE: Create a one dimensional gaussian kernel.
	 * NAME: Mike Heath
	 * DATE: 2/15/96
	 *******************************************************************************/
	void make_gaussian_kernel(float sigma, float *kernel, int *windowsize) {
		int i, center;
		float x, fx, sum = 0.0;

		*windowsize = 1 + 2 * ceil(2.5 * sigma);
		center = (*windowsize) / 2;

		if (VERBOSE)
			printf("      The kernel has %d elements.\n", *windowsize);

		for (i = 0; i < (*windowsize); i++) {
			x = (float) (i - center);
			fx = pow(2.71828, -0.5 * x * x / (sigma * sigma))
					/ (sigma * sqrt(6.2831853));
			kernel[i] = fx;
			sum += fx;
		}

		for (i = 0; i < (*windowsize); i++)
			kernel[i] /= sum;

		if (VERBOSE) {
			printf("The filter coefficients are:\n");
			for (i = 0; i < (*windowsize); i++)
				printf("kernel[%d] = %f\n", i, kernel[i]);
		}
	}

	void main(void) {
		// cout << "******In Make_GAUSSIAN_KERNEL**********" << endl;
		while (1) {
			wait(0, SC_MS);
			imagein.read(imgIn);

			make_gaussian_kernel(SIGMA, ker, &winsize);

			kernel.write(ker);
			windowsize.write(winsize);
			imageOut.write(imgIn);
		}

	}

	SC_CTOR(MAKE_GAUSSIAN_KERNEL) {
		winsize = 0;
		SC_THREAD(main);
		SET_STACK_SIZE
	}

};

SC_MODULE(BLUR_X) {

	sc_event e1, e2, e3, e4, data_received;

	IMAGE imgIn;
	KERNEL ker;
	TEMPIM tempim;
	int winsize;

	sc_fifo_in<IMAGE> image;
	sc_fifo_in<KERNEL> kernel;
	sc_fifo_in<int> windowsize;

	sc_fifo_out<TEMPIM> tpim;
	sc_fifo_out<KERNEL> kernel1;
	sc_fifo_out<int> windowsize1;

	void blurx1() {

		while (true) {
			wait(data_received);

			int r, c, cc; /* Counter variables. */
			float dot, sum;
			int rows = 1520;
			int cols = 2704;
			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the x - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the X-direction.\n");
			for (r = 0; r <= (rows / 4) * 1 - 1; r++) {
				for (c = 0; c < cols; c++) {
					dot = 0.0;
					sum = 0.0;
					for (cc = (-center); cc <= center; cc++) {
						if (((c + cc) >= 0) && ((c + cc) < cols)) {
							dot += (float) imgIn[r * cols + (c + cc)]
									* ker[center + cc];
							sum += ker[center + cc];
						}
					}
					tempim[r * cols + c] = dot / sum;
				}
			}
			wait(291 / 4, SC_MS);
			e1.notify(SC_ZERO_TIME);
		}
	}

	void blurx2() {

		while (true) {
			wait(data_received);

			int r, c, cc; /* Counter variables. */
			float dot, sum;
			int rows = 1520;
			int cols = 2704;
			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the x - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the X-direction.\n");
			for (r = (rows / 4) * 1; r <= (rows / 4) * 2 - 1; r++) {
				for (c = 0; c < cols; c++) {
					dot = 0.0;
					sum = 0.0;
					for (cc = (-center); cc <= center; cc++) {
						if (((c + cc) >= 0) && ((c + cc) < cols)) {
							dot += (float) imgIn[r * cols + (c + cc)]
									* ker[center + cc];
							sum += ker[center + cc];
						}
					}
					tempim[r * cols + c] = dot / sum;
				}
			}
			wait(291 / 4, SC_MS);
			e2.notify(SC_ZERO_TIME);
		}
	}

	void blurx3() {

		while (true) {

			wait(data_received);

			int r, c, cc; /* Counter variables. */
			float dot, sum;
			int rows = 1520;
			int cols = 2704;
			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the x - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the X-direction.\n");
			for (r = (rows / 4) * 2; r <= (rows / 4) * 3 - 1; r++) {
				for (c = 0; c < cols; c++) {
					dot = 0.0;
					sum = 0.0;
					for (cc = (-center); cc <= center; cc++) {
						if (((c + cc) >= 0) && ((c + cc) < cols)) {
							dot += (float) imgIn[r * cols + (c + cc)]
									* ker[center + cc];
							sum += ker[center + cc];
						}
					}
					tempim[r * cols + c] = dot / sum;
				}
			}
			wait(291 / 4, SC_MS);
			e3.notify(SC_ZERO_TIME);
		}
	}

	void blurx4() {

		while (true) {

			wait(data_received);

			int r, c, cc; /* Counter variables. */
			float dot, sum;
			int rows = 1520;
			int cols = 2704;
			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the x - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the X-direction.\n");
			for (r = (rows / 4) * 3; r <= (rows / 4) * 4 - 1; r++) {
				for (c = 0; c < cols; c++) {
					dot = 0.0;
					sum = 0.0;
					for (cc = (-center); cc <= center; cc++) {
						if (((c + cc) >= 0) && ((c + cc) < cols)) {
							dot += (float) imgIn[r * cols + (c + cc)]
									* ker[center + cc];
							sum += ker[center + cc];
						}
					}
					tempim[r * cols + c] = dot / sum;
				}
			}
			wait(291 / 4, SC_MS);
			e4.notify(SC_ZERO_TIME);
		}
	}

	void main(void) {
		// cout << "******In BLUR_X**********" << endl;
		while (1) {
			windowsize.read(winsize);
			kernel.read(ker);
			image.read(imgIn);

			data_received.notify(SC_ZERO_TIME);
			wait(e1 & e2 & e3 & e4);

			kernel1.write(ker);
			windowsize1.write(winsize);
			tpim.write(tempim);
		}
	}

	SC_CTOR(BLUR_X) {
		winsize = 0;
		SC_THREAD(main);
		SET_STACK_SIZE

		SC_THREAD(blurx1);
		SET_STACK_SIZE
		SC_THREAD(blurx2);
		SET_STACK_SIZE
		SC_THREAD(blurx3);
		SET_STACK_SIZE
		SC_THREAD(blurx4);
		SET_STACK_SIZE
	}

};

SC_MODULE(BLUR_Y) {

	sc_event e1, e2, e3, e4, data_received;

	KERNEL ker;
	TEMPIM tempim;
	SIMAGE smoothedim;
	int winsize;

	sc_fifo_in<KERNEL> kernel;
	sc_fifo_in<int> windowsize;
	sc_fifo_in<TEMPIM> tempIm;

	sc_fifo_out<SIMAGE> smoothedIm;

	void blury1() {

		while (true) {

			wait(data_received);

			int r, c, rr; /* Counter variables. */
			float dot, sum;

			int rows = 1520;
			int cols = 2704;

			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the y - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the Y-direction.\n");
			for (c = 0; c <= (cols / 4) * 1 - 1; c++) {
				for (r = 0; r < rows; r++) {
					sum = 0.0;
					dot = 0.0;
					for (rr = (-center); rr <= center; rr++) {
						if (((r + rr) >= 0) && ((r + rr) < rows)) {
							dot += tempim[(r + rr) * cols + c]
									* ker[center + rr];
							sum += ker[center + rr];
						}
					}
					smoothedim[r * cols + c] = (short int) (dot
							* BOOSTBLURFACTOR / sum + 0.5);
				}
			}
			wait(473 / 4, SC_MS);
			e1.notify(SC_ZERO_TIME);
		}
	}

	void blury2() {

		while (true) {

			wait(data_received);

			int r, c, rr; /* Counter variables. */
			float dot, sum;

			int rows = 1520;
			int cols = 2704;

			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the y - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the Y-direction.\n");
			for (c = (cols / 4) * 1; c <= (cols / 4) * 2 - 1; c++) {
				for (r = 0; r < rows; r++) {
					sum = 0.0;
					dot = 0.0;
					for (rr = (-center); rr <= center; rr++) {
						if (((r + rr) >= 0) && ((r + rr) < rows)) {
							dot += tempim[(r + rr) * cols + c]
									* ker[center + rr];
							sum += ker[center + rr];
						}
					}
					smoothedim[r * cols + c] = (short int) (dot
							* BOOSTBLURFACTOR / sum + 0.5);
				}
			}
			wait(473 / 4, SC_MS);
			e2.notify(SC_ZERO_TIME);
		}
	}

	void blury3() {

		while (true) {

			wait(data_received);

			int r, c, rr; /* Counter variables. */
			float dot, sum;

			int rows = 1520;
			int cols = 2704;

			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the y - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the Y-direction.\n");
			for (c = (cols / 4) * 2; c <= (cols / 4) * 3 - 1; c++) {
				for (r = 0; r < rows; r++) {
					sum = 0.0;
					dot = 0.0;
					for (rr = (-center); rr <= center; rr++) {
						if (((r + rr) >= 0) && ((r + rr) < rows)) {
							dot += tempim[(r + rr) * cols + c]
									* ker[center + rr];
							sum += ker[center + rr];
						}
					}
					smoothedim[r * cols + c] = (short int) (dot
							* BOOSTBLURFACTOR / sum + 0.5);
				}
			}
			wait(473 / 4, SC_MS);
			e3.notify(SC_ZERO_TIME);
		}
	}

	void blury4() {

		while (true) {

			wait(data_received);

			int r, c, rr; /* Counter variables. */
			float dot, sum;

			int rows = 1520;
			int cols = 2704;

			int center = winsize / 2;
			/****************************************************************************
			 * Blur in the y - direction.
			 ****************************************************************************/
			if (VERBOSE)
				printf("   Bluring the image in the Y-direction.\n");
			for (c = (cols / 4) * 3; c < (cols / 4) * 4 - 1; c++) {
				for (r = 0; r < rows; r++) {
					sum = 0.0;
					dot = 0.0;
					for (rr = (-center); rr <= center; rr++) {
						if (((r + rr) >= 0) && ((r + rr) < rows)) {
							dot += tempim[(r + rr) * cols + c]
									* ker[center + rr];
							sum += ker[center + rr];
						}
					}
					smoothedim[r * cols + c] = (short int) (dot
							* BOOSTBLURFACTOR / sum + 0.5);
				}
			}
			wait(473 / 4, SC_MS);
			e4.notify(SC_ZERO_TIME);
		}
	}

	void main(void) {
		// cout << "******In BLUR_Y**********" << endl;
		while (1) {
			windowsize.read(winsize);
			kernel.read(ker);
			tempIm.read(tempim);

			data_received.notify(SC_ZERO_TIME);
			wait(e1 & e2 & e3 & e4);

			smoothedIm.write(smoothedim);
		}
	}

	SC_CTOR(BLUR_Y) {
		winsize = 0;
		SC_THREAD(main);
		SET_STACK_SIZE

		SC_THREAD(blury1);
		SET_STACK_SIZE
		SC_THREAD(blury2);
		SET_STACK_SIZE
		SC_THREAD(blury3);
		SET_STACK_SIZE
		SC_THREAD(blury4);
		SET_STACK_SIZE
	}

};

SC_MODULE(GAUSSIAN_SMOOTH) {
	sc_fifo_in<IMAGE> pIn;
	sc_fifo_out<SIMAGE> pOut;

	sc_fifo<KERNEL> kernel;
	sc_fifo<KERNEL> kernel1;
	sc_fifo<int> windowsize;
	sc_fifo<int> windowsize1;
	sc_fifo<IMAGE> imageBuffer;
	sc_fifo<TEMPIM> tempim;

	MAKE_GAUSSIAN_KERNEL make_gaussian_kernel;
	BLUR_X blur_x;
	BLUR_Y blur_y;

	void before_end_of_elaboration() {
		make_gaussian_kernel.kernel.bind(kernel);
		make_gaussian_kernel.windowsize.bind(windowsize);
		make_gaussian_kernel.imagein.bind(pIn);
		make_gaussian_kernel.imageOut.bind(imageBuffer);

		blur_x.kernel.bind(kernel);
		blur_x.kernel1.bind(kernel1);
		blur_x.windowsize1.bind(windowsize1);
		blur_x.windowsize.bind(windowsize);
		blur_x.image.bind(imageBuffer);
		blur_x.tpim.bind(tempim);

		blur_y.kernel.bind(kernel1);
		blur_y.windowsize.bind(windowsize1);
		blur_y.tempIm.bind(tempim);
		blur_y.smoothedIm.bind(pOut);
	}

	SC_CTOR(GAUSSIAN_SMOOTH) :
			kernel("kernel", 1), kernel1("kernel1", 1), windowsize("windowsize",
					1), windowsize1("windowsize1", 1), imageBuffer(
					"imageBuffer", 1), tempim("tempim", 1), make_gaussian_kernel(
					"make_gaussian_kernel"), blur_x("blur_x"), blur_y("blur_y") {
		// cout << "****** In GAUSSIAN_SMOOTH **********" << endl;
	}

};

SC_MODULE(DERRIVATIVE_X_Y) {

	SIMAGE imageIn;

	SIMAGE imageOutX;
	SIMAGE imageOutY;

	sc_fifo_in<SIMAGE> pIn;
	sc_fifo_out<SIMAGE> deltaX;
	sc_fifo_out<SIMAGE> deltaY;

	void derrivative_x_y(short int *smoothedim, int rows, int cols,
			short int *delta_x, short int *delta_y) {
		int r, c, pos;

		if (VERBOSE)
			printf("   Computing the X-direction derivative.\n");
		for (r = 0; r < rows; r++) {
			pos = r * cols;
			delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
			pos++;
			for (c = 1; c < (cols - 1); c++, pos++) {
				delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
			}
			delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
		}

		if (VERBOSE)
			printf("   Computing the Y-direction derivative.\n");
		for (c = 0; c < cols; c++) {
			pos = c;
			delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
			pos += cols;
			for (r = 1; r < (rows - 1); r++, pos += cols) {
				delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
			}
			delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
		}
	}

	void main(void) {
		// cout << "****** In DERRIVATIVE_X_Y **********" << endl;
		while (1) {
			pIn.read(imageIn);
			wait(152, SC_MS);

			derrivative_x_y(imageIn, ROWS, COLS, imageOutX, imageOutY);

			deltaX.write(imageOutX);
			deltaY.write(imageOutY);
		}
	}

	SC_CTOR(DERRIVATIVE_X_Y) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}

};

SC_MODULE(MAGNITUDE_X_Y) {
	sc_fifo_in<SIMAGE> p1;
	sc_fifo_in<SIMAGE> p2;

	sc_fifo_out<SIMAGE> x;
	sc_fifo_out<SIMAGE> y;

	sc_fifo_out<SIMAGE> pOut;

	SIMAGE imageInX;
	SIMAGE imageInY;

	SIMAGE imageMag;

	void magnitude_x_y(short int *delta_x, short int *delta_y, int rows,
			int cols, short int *magnitude) {
		int r, c, pos, sq1, sq2;

		for (r = 0, pos = 0; r < rows; r++) {
			for (c = 0; c < cols; c++, pos++) {
				sq1 = (int) delta_x[pos] * (int) delta_x[pos];
				sq2 = (int) delta_y[pos] * (int) delta_y[pos];
				magnitude[pos] =
						(short) (0.5 + sqrt((float) sq1 + (float) sq2));
			}
		}

	}

	void main(void) {
		// cout << "****** In MAGNITUDE_X_Y **********" << endl;
		while (1) {
			p1.read(imageInX);
			p2.read(imageInY);
			wait(156, SC_MS);

			x.write(imageInX);
			y.write(imageInY);

			magnitude_x_y(imageInX, imageInY, ROWS, COLS, imageMag);

			pOut.write(imageMag);
		}
	}

	SC_CTOR(MAGNITUDE_X_Y) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}

};

SC_MODULE(NON_MAX_SUPPLY) {

	sc_fifo_in<SIMAGE> pIn;
	sc_fifo_in<SIMAGE> dx;
	sc_fifo_in<SIMAGE> dy;

	sc_fifo_out<IMAGE> pOut;
	sc_fifo_out<SIMAGE> pOut2;

	SIMAGE imageX;
	SIMAGE imageY;
	SIMAGE imageMag;

	IMAGE imageNMS;

	void non_max_supp(short *mag, short *gradx, short *grady, int nrows,
			int ncols, unsigned char *result) {
		int rowcount, colcount, count;
		short *magrowptr, *magptr;
		short *gxrowptr, *gxptr;
		short *gyrowptr, *gyptr, z1, z2;
		short m00, gx, gy;
		float mag1, mag2, xperp, yperp;
		unsigned char *resultrowptr, *resultptr;

		for (count = 0, resultrowptr = result, resultptr = result
				+ ncols * (nrows - 1); count < ncols;
				resultptr++, resultrowptr++, count++) {
			*resultrowptr = *resultptr = (unsigned char) 0;
		}

		for (count = 0, resultptr = result, resultrowptr = result + ncols - 1;
				count < nrows; count++, resultptr += ncols, resultrowptr +=
						ncols) {
			*resultptr = *resultrowptr = (unsigned char) 0;
		}

		for (rowcount = 1, magrowptr = mag + ncols + 1, gxrowptr = gradx + ncols
				+ 1, gyrowptr = grady + ncols + 1, resultrowptr = result + ncols
				+ 1; rowcount <= nrows - 2;	// bug fix 3/29/17, RD
				rowcount++, magrowptr += ncols, gyrowptr += ncols, gxrowptr +=
						ncols, resultrowptr += ncols) {
			for (colcount = 1, magptr = magrowptr, gxptr = gxrowptr, gyptr =
					gyrowptr, resultptr = resultrowptr; colcount <= ncols - 2;	// bug fix 3/29/17, RD
					colcount++, magptr++, gxptr++, gyptr++, resultptr++) {
				m00 = *magptr;
				if (m00 == 0) {
					*resultptr = (unsigned char) NOEDGE;
				} else {
					xperp = -(gx = *gxptr) / ((float) m00);
					yperp = (gy = *gyptr) / ((float) m00);
					// gx = *gxptr;
					// gy = *gyptr;
					// xperp = -(gx << 16) / m00;
					// yperp = (gy << 16) / m00;
				}

				if (gx >= 0) {
					if (gy >= 0) {
						if (gx >= gy) {
							/* 111 */
							/* Left point */
							z1 = *(magptr - 1);
							z2 = *(magptr - ncols - 1);

							mag1 = (m00 - z1) * xperp + (z2 - z1) * yperp;

							/* Right point */
							z1 = *(magptr + 1);
							z2 = *(magptr + ncols + 1);

							mag2 = (m00 - z1) * xperp + (z2 - z1) * yperp;
						} else {
							/* 110 */
							/* Left point */
							z1 = *(magptr - ncols);
							z2 = *(magptr - ncols - 1);

							mag1 = (z1 - z2) * xperp + (z1 - m00) * yperp;

							/* Right point */
							z1 = *(magptr + ncols);
							z2 = *(magptr + ncols + 1);

							mag2 = (z1 - z2) * xperp + (z1 - m00) * yperp;
						}
					} else {
						if (gx >= -gy) {
							/* 101 */
							/* Left point */
							z1 = *(magptr - 1);
							z2 = *(magptr + ncols - 1);

							mag1 = (m00 - z1) * xperp + (z1 - z2) * yperp;

							/* Right point */
							z1 = *(magptr + 1);
							z2 = *(magptr - ncols + 1);

							mag2 = (m00 - z1) * xperp + (z1 - z2) * yperp;
						} else {
							/* 100 */
							/* Left point */
							z1 = *(magptr + ncols);
							z2 = *(magptr + ncols - 1);

							mag1 = (z1 - z2) * xperp + (m00 - z1) * yperp;

							/* Right point */
							z1 = *(magptr - ncols);
							z2 = *(magptr - ncols + 1);

							mag2 = (z1 - z2) * xperp + (m00 - z1) * yperp;
						}
					}
				} else {
					if ((gy = *gyptr) >= 0) {
						if (-gx >= gy) {
							/* 011 */
							/* Left point */
							z1 = *(magptr + 1);
							z2 = *(magptr - ncols + 1);

							mag1 = (z1 - m00) * xperp + (z2 - z1) * yperp;

							/* Right point */
							z1 = *(magptr - 1);
							z2 = *(magptr + ncols - 1);

							mag2 = (z1 - m00) * xperp + (z2 - z1) * yperp;
						} else {
							/* 010 */
							/* Left point */
							z1 = *(magptr - ncols);
							z2 = *(magptr - ncols + 1);

							mag1 = (z2 - z1) * xperp + (z1 - m00) * yperp;

							/* Right point */
							z1 = *(magptr + ncols);
							z2 = *(magptr + ncols - 1);

							mag2 = (z2 - z1) * xperp + (z1 - m00) * yperp;
						}
					} else {
						if (-gx > -gy) {
							/* 001 */
							/* Left point */
							z1 = *(magptr + 1);
							z2 = *(magptr + ncols + 1);

							mag1 = (z1 - m00) * xperp + (z1 - z2) * yperp;

							/* Right point */
							z1 = *(magptr - 1);
							z2 = *(magptr - ncols - 1);

							mag2 = (z1 - m00) * xperp + (z1 - z2) * yperp;
						} else {
							/* 000 */
							/* Left point */
							z1 = *(magptr + ncols);
							z2 = *(magptr + ncols + 1);

							mag1 = (z2 - z1) * xperp + (m00 - z1) * yperp;

							/* Right point */
							z1 = *(magptr - ncols);
							z2 = *(magptr - ncols - 1);

							mag2 = (z2 - z1) * xperp + (m00 - z1) * yperp;
						}
					}
				}

				if ((mag1 > 0.0) || (mag2 > 0.0)) {
					*resultptr = (unsigned char) NOEDGE;
				} else {
					if (mag2 == 0.0)
						*resultptr = (unsigned char) NOEDGE;
					else
						*resultptr = (unsigned char) POSSIBLE_EDGE;
				}
			}
		}
	}

	void main(void) {
		// cout << "****** In NON_MAX_SUPPLY **********" << endl;
		while (1) {
			pIn.read(imageMag);
			dx.read(imageX);
			dy.read(imageY);
			wait(638, SC_MS);

			non_max_supp(imageMag, imageX, imageY, ROWS, COLS, imageNMS);
			pOut.write(imageNMS);
			pOut2.write(imageMag);
		}

	}

	SC_CTOR(NON_MAX_SUPPLY) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}

};

SC_MODULE(APPLY_HYSTERESIS) {

	sc_fifo_in<SIMAGE> p1;
	sc_fifo_in<IMAGE> p2;

	sc_fifo_out<IMAGE> pOut;

	SIMAGE imageMAG;
	IMAGE imageNMS;

	IMAGE imageOUT;

	void follow_edges(unsigned char *edgemapptr, short *edgemagptr,
			short lowval, int cols) {
		short *tempmagptr;
		unsigned char *tempmapptr;
		int i;
		int x[8] = { 1, 1, 0, -1, -1, -1, 0, 1 }, y[8] = { 0, 1, 1, 1, 0, -1,
				-1, -1 };

		for (i = 0; i < 8; i++) {
			tempmapptr = edgemapptr - y[i] * cols + x[i];
			tempmagptr = edgemagptr - y[i] * cols + x[i];

			if ((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)) {
				*tempmapptr = (unsigned char) EDGE;
				follow_edges(tempmapptr, tempmagptr, lowval, cols);
			}
		}
	}

	void apply_hysteresis(short int *mag, unsigned char *nms, int rows,
			int cols, float tlow, float thigh, unsigned char *edge) {
		int r, c, pos, numedges, highcount, lowthreshold, highthreshold,
				hist[32768];
		short int maximum_mag;

		for (r = 0, pos = 0; r < rows; r++) {
			for (c = 0; c < cols; c++, pos++) {
				if (nms[pos] == POSSIBLE_EDGE)
					edge[pos] = POSSIBLE_EDGE;
				else
					edge[pos] = NOEDGE;
			}
		}

		for (r = 0, pos = 0; r < rows; r++, pos += cols) {
			edge[pos] = NOEDGE;
			edge[pos + cols - 1] = NOEDGE;
		}
		pos = (rows - 1) * cols;
		for (c = 0; c < cols; c++, pos++) {
			edge[c] = NOEDGE;
			edge[pos] = NOEDGE;
		}

		for (r = 0; r < 32768; r++)
			hist[r] = 0;
		for (r = 0, pos = 0; r < rows; r++) {
			for (c = 0; c < cols; c++, pos++) {
				if (edge[pos] == POSSIBLE_EDGE)
					hist[mag[pos]]++;
			}
		}

		for (r = 1, numedges = 0; r < 32768; r++) {
			if (hist[r] != 0)
				maximum_mag = r;
			numedges += hist[r];
		}

		highcount = (int) (numedges * thigh + 0.5);

		r = 1;
		numedges = hist[1];
		while ((r < (maximum_mag - 1)) && (numedges < highcount)) {
			r++;
			numedges += hist[r];
		}
		highthreshold = r;
		lowthreshold = (int) (highthreshold * tlow + 0.5);

		if (VERBOSE) {
			printf(
					"The input low and high fractions of %f and %f computed to\n",
					tlow, thigh);
			printf("magnitude of the gradient threshold values of: %d %d\n",
					lowthreshold, highthreshold);
		}

		for (r = 0, pos = 0; r < rows; r++) {
			for (c = 0; c < cols; c++, pos++) {
				if ((edge[pos] == POSSIBLE_EDGE)
						&& (mag[pos] >= highthreshold)) {
					edge[pos] = EDGE;
					follow_edges((edge + pos), (mag + pos), lowthreshold, cols);
				}
			}
		}

		for (r = 0, pos = 0; r < rows; r++) {
			for (c = 0; c < cols; c++, pos++)
				if (edge[pos] != EDGE)
					edge[pos] = NOEDGE;
		}
	}

	void main(void) {
		while (1) {
			// cout << "****** In APPLY_HYSTERESIS **********" << endl;
			p1.read(imageMAG);
			p2.read(imageNMS);
			wait(219, SC_MS);

			apply_hysteresis(imageMAG, imageNMS, ROWS, COLS, TLOW, THIGH,
					imageOUT);

			pOut.write(imageOUT);
		}
	}
	SC_CTOR(APPLY_HYSTERESIS) {
		SC_THREAD(main);
		SET_STACK_SIZE
	}

};

SC_MODULE(DUT) {

	sc_fifo_in<IMAGE> p1;
	sc_fifo_out<IMAGE> p2;

	sc_fifo<SIMAGE> smoothdim;
	sc_fifo<SIMAGE> deltaX;
	sc_fifo<SIMAGE> deltaY;
	sc_fifo<SIMAGE> dX;
	sc_fifo<SIMAGE> dY;
	sc_fifo<SIMAGE> magnitude;
	sc_fifo<SIMAGE> magnitude2;

	sc_fifo<IMAGE> nms;

	GAUSSIAN_SMOOTH gaussian_smooth;
	DERRIVATIVE_X_Y derrivative_x_y;
	MAGNITUDE_X_Y magnitude_x_y;
	NON_MAX_SUPPLY non_max_supply;
	APPLY_HYSTERESIS apply_hysteresis;

	void before_end_of_elaboration() {
		gaussian_smooth.pIn.bind(p1);
		gaussian_smooth.pOut.bind(smoothdim);

		derrivative_x_y.pIn.bind(smoothdim);
		derrivative_x_y.deltaX.bind(deltaX);
		derrivative_x_y.deltaY.bind(deltaY);

		magnitude_x_y.p1.bind(deltaX);
		magnitude_x_y.p2.bind(deltaY);
		magnitude_x_y.x.bind(dX);
		magnitude_x_y.y.bind(dY);
		magnitude_x_y.pOut.bind(magnitude);

		non_max_supply.pIn.bind(magnitude);
		non_max_supply.dx.bind(dX);
		non_max_supply.dy.bind(dY);
		non_max_supply.pOut.bind(nms);
		non_max_supply.pOut2.bind(magnitude2);

		apply_hysteresis.p1.bind(magnitude2);
		apply_hysteresis.p2.bind(nms);
		apply_hysteresis.pOut.bind(p2);

	}

	SC_CTOR(DUT) :
			smoothdim("smoothdim", 1), deltaX("deltaX", 1), deltaY("deltaY", 1), dX(
					"dX", 1), dY("dY", 1), magnitude("magnitude", 1), magnitude2(
					"magnitude2", 1), nms("nms", 1), gaussian_smooth(
					"gaussian_smooth"), derrivative_x_y("derrivative_x_y"), magnitude_x_y(
					"magnitude_x_y"), non_max_supply("non_max_supply"), apply_hysteresis(
					"apply_hysteresis") {

		// cout << "******In DUT**********" << endl;

	}
};

SC_MODULE(Platform) {
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;

	DataIn din;
	DUT canny;
	DataOut dout;

	void before_end_of_elaboration() {
		din.ImgIn.bind(ImgIn);
		din.ImgOut.bind(q1);
		canny.p1.bind(q1);
		canny.p2.bind(q2);
		dout.ImgIn.bind(q2);
		dout.ImgOut.bind(ImgOut);
	}

	SC_CTOR(Platform) :
			q1("q1", 1), q2("q2", 1), din("din"), canny("canny"), dout("dout") {
		// cout << "****** In Platform **********" << endl;
	}
};

SC_MODULE(Top) {
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;
	sc_fifo<sc_time> frameTimeStamp;
	Stimulus stimulus;
	Platform platform;
	Monitor monitor;

	void before_end_of_elaboration() {
		stimulus.ImgOut.bind(q1);
		stimulus.frameTimeOut.bind(frameTimeStamp);
		platform.ImgIn.bind(q1);
		platform.ImgOut.bind(q2);
		monitor.ImgIn.bind(q2);
		monitor.frameTimeIn.bind(frameTimeStamp);
	}

	SC_CTOR(Top) :
			q1("q1", 1), q2("q2", 1), frameTimeStamp("frameTimeStamp",
			AVAIL_IMG), stimulus("stimulus"), platform("platform"), monitor(
					"monitor") {
		//	cout << "****** In Top **********" << endl;
	}
};

Top top("top");

int sc_main(int argc, char* argv[]) {
	// cout << "****** In sc_main **********" << endl;
	sc_start();
	return 0;
}
