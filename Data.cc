//////////////////////////////////////////////////////////
//      Rémi Cura , Thales IGN, 08/2014                 //
//          Street Gen snapping                         //
//                                                      //
//////////////////////////////////////////////////////////
/**
  This file holds the Data.h class implementation

  */

#include "Data.h";

/** Constructor for data input/storage
  @param name of the file containing the input data : the format is very specific
  @param name of the file where to write the temporary result for each iteration
  */
DataStorage::DataStorage(const std::string& i_filename, const std::string& o_filename ) {
    readData(i_filename);
    }


DataStorage::readData(const std::string& i_filename){
    //opening files
    FILE* i_fptr = fopen(i_filename.c_str(), "r");

    if (i_fptr == NULL) {
      stder << "Error: unable to open file " << i_filename;
      exit(1);
    };

    //reading data parameter : how much we have to read after this.
        // see header of this class to get an idea of the required input file format
        FscanfOrDie(i_fptr, "%d", &num_nodes_);
        FscanfOrDie(i_fptr, "%d", &num_edges_);
        FscanfOrDie(i_fptr, "%d", &num_observations_);

        std::cout << "num nodes : " << num_nodes() << "\n"
                  << "num edges : " << num_edges() << "\n"
                  << "num edges : " << num_observations() << "\n";

    //reading data
        char line[1000];

    //reading the nodes data :


    //node_id::int;X::double;Yi_filename::double;Z::double;is_in_intersection::int
    for (int i = 0; i < num_nodes_; ++i) {
        fgets(line, sizeof line, i_fptr);
        node t_node = new node();
        double coor = new double[3];
        if (*line == '#') continue; /* ignore comment line */
        if (sscanf(line, "%+d;%G;%G;%G;%+d",t_node.node_id,coor[0],coor[1],coor[2],t_node.is_in_intersection) != 5) {
            std::cerr<< "error when trying ot read the nodes, wrong format" ;
            exit(1);
        } else {
            t_node.position = coor;
            nodes_[i] = t_node;
        }
    }

    //reading the edge data :




    // This wil die horribly on invalid files. Them's the breaks.
    FscanfOrDie(i_fptr, "%d", &num_cameras_);
    FscanfOrDie(i_fptr, "%d", &num_points_);
    FscanfOrDie(i_fptr, "%d", &num_observations_);

    VLOG(1) << "Header: " << num_cameras_
            << " " << num_points_
            << " " << num_observations_;

    point_index_ = new int[num_observations_];
    camera_index_ = new int[num_observations_];
    observations_ = new double[2 * num_observations_];

    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
    parameters_ = new double[num_parameters_];

    for (int i = 0; i < num_observations_; ++i) {
      FscanfOrDie(fptr, "%d", camera_index_ + i);
      FscanfOrDie(fptr, "%d", point_index_ + i);
      for (int j = 0; j < 2; ++j) {
        FscanfOrDie(fptr, "%lf", observations_ + 2*i + j);
      }
    }

    for (int i = 0; i < num_parameters_; ++i) {
      FscanfOrDie(fptr, "%lf", parameters_ + i);
    }

    fclose(fptr);

    use_quaternions_ = use_quaternions;
    if (use_quaternions) {
      // Switch the angle-axis rotations to quaternions.
      num_parameters_ = 10 * num_cameras_ + 3 * num_points_;
      double* quaternion_parameters = new double[num_parameters_];
      double* original_cursor = parameters_;
      double* quaternion_cursor = quaternion_parameters;
      for (int i = 0; i < num_cameras_; ++i) {
        AngleAxisToQuaternion(original_cursor, quaternion_cursor);
        quaternion_cursor += 4;
        original_cursor += 3;
        for (int j = 4; j < 10; ++j) {
         *quaternion_cursor++ = *original_cursor++;
        }
      }
      // Copy the rest of the points.
      for (int i = 0; i < 3 * num_points_; ++i) {
        *quaternion_cursor++ = *original_cursor++;
      }
      // Swap in the quaternion parameters.
      delete []parameters_;
      parameters_ = quaternion_parameters;

  }

}
