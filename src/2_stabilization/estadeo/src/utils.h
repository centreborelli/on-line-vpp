#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> 


size_t read_video(
  char *name,       //file name
  unsigned char *I, //video to read
  int size          //size of the video
);

size_t write_video(
  char *name,       //file name
  unsigned char *I, //video to write
  int size          //size of the video
);

void save_transform(
  char  *name,   //file name
  float *H,      //transformation
  int   nparams  //number of parameters of the transformations
);


//class for recording runtimes of the method
class Timer {

  public:
    Timer():avg1(), avg2(), avg3() {}
    void set_t1(){gettimeofday(&t1, NULL);}
    void set_t2(){gettimeofday(&t2, NULL);}
    void set_t3(){gettimeofday(&t3, NULL);}
    void set_t4(){gettimeofday(&t4, NULL); avg_update();}
    
    void avg_update(){
      avg1+=((t2.tv_sec-t1.tv_sec)*1000000u+t2.tv_usec-t1.tv_usec)/1.e6;
      avg2+=((t3.tv_sec-t2.tv_sec)*1000000u+t3.tv_usec-t2.tv_usec)/1.e6;
      avg3+=((t4.tv_sec-t3.tv_sec)*1000000u+t4.tv_usec-t3.tv_usec)/1.e6;
    }
    
    void print_time(int f) {
      printf(" Processing frame %d: T(%.4fs, %.7fs, %.4fs) \n", f, 
            ((t2.tv_sec-t1.tv_sec)*1000000u+t2.tv_usec-t1.tv_usec)/1.e6,
            ((t3.tv_sec-t2.tv_sec)*1000000u+t3.tv_usec-t2.tv_usec)/1.e6,
            ((t4.tv_sec-t3.tv_sec)*1000000u+t4.tv_usec-t3.tv_usec)/1.e6);
    }
    
    void print_avg_time(int nframes) {
      float total=avg1+avg2+avg3;
      float average=(avg1+avg2+avg3)/nframes;
      printf(
        "\n Average time per frame: %.4fs -> T(%.4fs, %.7fs, %.4fs) \n", 
        average, avg1/nframes, avg2/nframes, avg3/nframes
      ); 
      printf(
        "\n Total time: %.4fs -> T(%.4fs, %.7fs, %.4fs) \n", 
        total, avg1, avg2, avg3
      ); 
    }
    
  private:
    float avg1, avg2, avg3;
    struct timeval t1, t2, t3, t4;  
};


#endif
