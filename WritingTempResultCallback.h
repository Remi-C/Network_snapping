#ifndef WRITING_TEMP_RESULT_CALLBACK_H_
#define WRITING_TEMP_RESULT_CALLBACK_H_
//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////

/** This class hold the functor that will output the intermediary result at each iteration.

  we write in a text file, using the csv format,
  1 line per edge per iteration
  1 id per edge per iteration
  we use WKT for geometry
      id;geom;start;end;
      1;LINESTRING( 0 0 0, 6 0 0) ;2014-08-29 13:00:00;2014-08-29 13:00:50
      1;LINESTRING( 0.5 0.5 0, 6 0.5 0) ;2014-08-29 13:01:00;2014-08-29 13:01:50
      1;LINESTRING( 1 1 0, 6 1 0) ;2014-08-29 13:02:00;2014-08-29 13:02:50
  */

#include <string>  // we need to store the file name
#include <fstream>  // we need some function to write in files
using namespace std; //maybe we could limit the visibility?

class WritingTempResultCallback : public ceres::IterationCallback {
 public:

  //! constructor : open the file, write the csv header in it
  WritingTempResultCallback(string file_n,int i_)
      : file_name(file_n) ,  i(i_) {
      file.open(file_n.c_str()) ;
      file << "id;geom;start;end"<<"\n";
  }

  ~WritingTempResultCallback(){
      //closing the connection to file system :
      file.close();
  }

  virtual ceres::CallbackReturnType operator(
    )(
          const ceres::IterationSummary& summary
     //arg for the functor. none needed.
          ) {
    //writing edge with updated node position in the file :
      ++i;
      file << " i : " << i << "some values : "<<"\n" ;
//      for (int k = 0; k <jNumNodes; ++k ){
//       std::cout << "NEw n_"<<k<<": \n"
//                   << node_position[3 * k] <<": \n"
//                  << node_position[3 * k+1]  <<": \n"
//                  << node_position[3 * k+2]
//                 << "\n ";
        return ceres::SOLVER_CONTINUE;
     }

 private:
   const string file_name;
   ofstream file;
   int i;

};


#endif //end of #define