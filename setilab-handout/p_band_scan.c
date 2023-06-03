#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <sched.h>
#include <unistd.h>
#include <pthread.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0


// TODO Create a struct that contains all of the necessary arguments to
// do the computation
typedef struct args {
  // Analyze_signal arguments
  signal* sig;
  int filter_order;
  int num_bands;
  double* lb;
  double* ub;
  int id;

  // Derived Quants in analyze_signal
  int bandwidth;
  double* band_power;

  // Main arguments
  int num_threads;
  int num_processors;

  int* num_bands_per_thread;
  int* final_array;

  double* filter_coefficients;


  // double band_power[];
} args;

args swag_1;

void usage() {
  printf("usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {
  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}


// TODO Create the worker function. Note that we are ~not~ going
// to be passing in any ~actual~ arguments into void. The entire
// point of this is because whenever we are calling the pthread_create()
// function, we are going to be passing in the struct which ontains
// all of the differnet fields that are necessary for the computation
// that we want the thread to actually be completing.
void* worker (void* arg)  {

  args *swag = (args*) arg; // Type-casting the struct so that we can utilize it
                          //
  //
  for (int band = swag_1.final_array[swag_1.id]; band < swag_1.final_array[swag_1.id + 1]; band++) {
    // Make the filter
    generate_band_pass(swag->sig->Fs,
                       band * swag->bandwidth + 0.0001, // keep within limits
                       (band + 1) * swag->bandwidth - 0.0001,
                       swag->filter_order,
                       swag->filter_coefficients); //
    hamming_window(swag->filter_order, swag->filter_coefficients);


    // Convolve
    convolve_and_compute_power(swag->sig->num_samples,
                               swag->sig->data,
                               swag->filter_order,
                               swag->filter_coefficients,
                               &(swag->band_power[band]));
  }
  pthread_exit(NULL);
}


// TODO
// - Create an array that initializes the thread s
// - Pass in the arguments "num_threads" and "num_processors"
int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub, int num_threads, int num_processors) {

  // Initialize the memory addresses for which the threads will
  // perform their computations
  pthread_t threads_ids[num_threads];

  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  double band_power[num_bands];

  // Declare the struct that contains all of the arguments
  // necessary for splitting all of the threads
  //

  // Actually initialize the threads and give them the space within the thread_ids array
  //
  // Creating a dynamic array that allocates memory for any number of processor threads
  args swag_2[num_threads];
  pthread_t* thread_array = malloc(sizeof(pthread_t) * num_threads);


  // NOTE: pthread_create(thread_handle, attirbutes, thread_function, function_argument)
  for (int i = 0; i < num_threads; i++) {
    // swag_2[i].sig = sig;
    swag_2[i].filter_order = filter_order;
    // swag_2[i].num_bands = num_bands;
    // swag_2[i].lb = lb;
    // swag_2[i].ub = ub;
    swag_2[i].bandwidth = bandwidth;
    swag_2[i].band_power = band_power;
    swag_2[i].num_threads = num_threads;
    // swag_2[i].num_processors = num_processors;
    swag_2[i].id = i;
    swag_2[i].filter_coefficients = (double*) malloc(sizeof(double) * (filter_order + 1));
    pthread_create( &(thread_array[i]), NULL, worker, (void*) &swag_2[i]);
  } 

  for (long i = 0; i < num_threads; i++) {
    pthread_join(threads_ids[i], NULL);
  }

  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n"
         "User time        %lf seconds\n"
         "System time      %lf seconds\n"
         "Page faults      %ld\n"
         "Page swaps       %ld\n"
         "Blocks of I/O    %ld\n"
         "Signals caught   %ld\n"
         "Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

// TODO Allow for the pass in of arguments
int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  int num_threads = atoi(argv[6]);
  int num_processors = atoi(argv[7]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n"
         "file:     %s\n"
         "Fs:       %lf Hz\n"
         "order:    %d\n"
         "bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;

  int var = num_bands / num_threads;
  int remainder = num_bands % num_threads;

  swag_1.num_bands_per_thread = (int*)malloc(sizeof(int) * num_threads);
  swag_1.final_array = (int*)malloc(sizeof(int) * (num_threads + 1));

  for (int i = 0; i < num_threads; i++) {
    swag_1.num_bands_per_thread[i] = var;
  }

  for (int j = 0; j < remainder; j++) {
    swag_1.num_bands_per_thread[j] += 1;
  }

  swag_1.final_array[0] = 0;
  for (int k = 1; k < num_threads + 1; k++) {
    swag_1.final_array[k] = swag_1.num_bands_per_thread[k-1] + swag_1.final_array[k-1];
  }

  // NOTE: Added the "num_threads" and the "num_processors" argument to the analyze_signal functoi n
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_threads, num_processors)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
