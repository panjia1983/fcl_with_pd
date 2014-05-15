#define BOOST_TEST_MODULE "FCL_XMLDATA"
#include <tinyxml.h>

#include <boost/test/unit_test.hpp>

#include "fcl/intersect.h"
#include "fcl/collision.h"
#include "fcl/math/sampling.h"
#include "fcl/BVH/BVH_model.h"
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <sstream>
#include <fstream>


#include "libsvm_classifier.h"
#include "fcl/penetration_depth.h"
#include "fcl/collision_data.h"
#include "fcl_resources/config.h"
#include "test_fcl_utility.h"


using namespace fcl;

static void loadSceneFile(const std::string& filename,
                          std::vector<std::vector<Vec3f> >& points_array,
                          std::vector<std::vector<Triangle> >& triangles_array,
                          std::vector<std::pair<Transform3f, Transform3f> >& motions)
{
  TiXmlDocument doc(filename.c_str());
  if(doc.LoadFile())
  {
    TiXmlHandle hdoc(&doc);
    TiXmlElement* model, *motion;
    model = doc.FirstChildElement("MODEL");
    if(model)
    {
      TiXmlElement* object = model->FirstChildElement("OBJECT");
      int i = 0;
      while(object)
      {
        int object_id = -1;
        object->Attribute("id", &object_id);
        // std::cout << object_id << std::endl;

        std::vector<Vec3f> points;
        std::vector<Triangle> triangles;

        TiXmlElement* grid = object->FirstChildElement("GRID");
        int n_vertices = 0;
        while(grid)
        {
          int grid_id;
          double grid_x, grid_y, grid_z;

          grid->Attribute("id", &grid_id);
          grid->Attribute("x", &grid_x);
          grid->Attribute("y", &grid_y);
          grid->Attribute("z", &grid_z);

          Vec3f p(grid_x, grid_y, grid_z);

          if(grid_id - 1 == (int)points.size())
            points.push_back(p);
          else if(grid_id - 1 < (int)points.size())
            points[grid_id - 1] = p;
          else // if(grid_id - 1 > points.size())
          {
            points.resize(grid_id);
            points.back() = p;
          }
                          
          n_vertices++;
          grid = grid->NextSiblingElement("GRID");
        }

        // std::cout << "#vertices " << n_vertices << std::endl;

        TiXmlElement* tri = object->FirstChildElement("TRIA");
        int n_tris = 0;
        while(tri)
        {
          int tri_id;
          int v1, v2, v3;

          tri->Attribute("id", &tri_id);
          tri->Attribute("g1", &v1);
          tri->Attribute("g2", &v2);
          tri->Attribute("g3", &v3);

          Triangle t(v1-1, v2-1, v3-1);

          if(tri_id - 1 == (int)triangles.size())
            triangles.push_back(t);
          else if(tri_id - 1 < (int)triangles.size())
            triangles[tri_id - 1] = t;
          else
          {
            triangles.resize(tri_id);
            triangles.back() = t;
          }
          
          n_tris++;
          tri = tri->NextSiblingElement("TRIA");
        }

        // std::cout << "#triangles " << n_tris << std::endl;

        if(object_id - 1 == (int)points_array.size())
        {
          points_array.push_back(points);
          triangles_array.push_back(triangles);
        }
        else if(object_id - 1 < (int)points_array.size())
        {
          points_array[object_id - 1] = points;
          triangles_array[object_id - 1] = triangles;
        }
        else
        {
          points_array.resize(object_id);
          triangles_array.resize(object_id);
          points_array.back() = points;
          triangles_array.back() = triangles;
        }

        object = object->NextSiblingElement("OBJECT");
        i++;
      }

      // std::cout << "#objects " << i << std::endl;
    }

    motion = doc.FirstChildElement("MOTION");
    if(motion)
    {
      TiXmlElement* frame = motion->FirstChildElement("FRAME");
      int n_frame = 0;
      while(frame)
      {
        int frame_id;
        double frame_time;
        frame->Attribute("id", &frame_id);
        frame->Attribute("time", &frame_time);


        Vec3f T1, T2;
        Matrix3f R1, R2;
        const char* obj1_pos_string = frame->Attribute("obj1_pos");
        const char* obj2_pos_string = frame->Attribute("obj2_pos");
        const char* obj1_dc_string = frame->Attribute("obj1_dc");
        const char* obj2_dc_string = frame->Attribute("obj2_dc");

        std::stringstream s1_pos(obj1_pos_string);
        s1_pos >> T1[0] >> T1[1] >> T1[2];
        std::stringstream s2_pos(obj2_pos_string);
        s2_pos >> T2[0] >> T2[1] >> T2[2];
        std::stringstream s1_mat(obj1_dc_string);
        for(int j = 0; j < 3; ++j)
          for(int k = 0; k < 3; ++k)
            s1_mat >> R1(j, k);
        std::stringstream s2_mat(obj2_dc_string);
        for(int j = 0; j < 3; ++j)
          for(int k = 0; k < 3; ++k)
            s2_mat >> R2(j, k);


        std::pair<Transform3f, Transform3f> motion = std::make_pair(Transform3f(R1, T1), Transform3f(R2, T2));
        if(frame_id - 1 == (int)motions.size())
          motions.push_back(motion);
        if(frame_id - 1 < (int)motions.size())
          motions[frame_id - 1] = motion;
        else
          motions.push_back(motion);        
        
        frame = frame->NextSiblingElement("FRAME");
        n_frame++;
      }

      // std::cout << "#frames " << n_frame << std::endl;
    }
  }
  else
    std::cerr << "Failed to load file " << filename << std::endl;    
}

BOOST_AUTO_TEST_CASE(scene_test)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;
  boost::filesystem::path path(TEST_RESOURCES_DIR);
  std::string filename = (path / "scenario-1-2-3/Model_1_Scenario_3.txt").string();
  loadSceneFile(filename, points_array, triangles_array, motions);

  BVHModel<OBBRSS> m1;
  BVHModel<OBBRSS> m2;
  m1.beginModel();
  m1.addSubModel(points_array[0], triangles_array[0]);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(points_array[1], triangles_array[1]);
  m2.endModel();

  
  CollisionResult result0;
  CollisionRequest request0;
  collide(&m1, motions[0].first, &m2, motions[0].second, request0, result0);
  // std::cout << result0.numContacts() << std::endl;

  CollisionResult result1;
  CollisionRequest request1;
  collide(&m1, motions[1].first, &m2, motions[1].second, request1, result1);
  // std::cout << result1.numContacts() << std::endl;
}

static std::string num2string(int i)
{
  std::string result;
  std::ostringstream convert;
  convert << i;
  result = convert.str();
  return result;
}

static void xml2obj(const std::string& in_filename, const std::string& out_filename_base)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;
  loadSceneFile(in_filename, points_array, triangles_array, motions);

  std::size_t n_obj = points_array.size();
  // save objs in local frame
  for(std::size_t i = 0; i < n_obj; ++i)
  {
    std::string out_filenameL = out_filename_base + num2string(i+1) + "L.obj";
    saveOBJFile(out_filenameL.c_str(), points_array[i], triangles_array[i]);
  }

  /*
  // save objs in frame 1
  for(std::size_t i = 0; i < n_obj; ++i)
  {
    std::string out_filenameF1 = out_filename_base + num2string(i+1) + "Frame1.obj";
    std::vector<Vec3f> points(points_array[i].size());
    for(std::size_t j = 0; j < points.size(); ++j)
      points[j] = motions[i].first.transform(points_array[i][j]);

    saveOBJFile(out_filenameF1.c_str(), points, triangles_array[i]);
  }

  // save objs in frame 2
  for(std::size_t i = 0; i < n_obj; ++i)
  {
    std::string out_filenameF2 = out_filename_base + num2string(i+1) + "Frame2.obj";
    std::vector<Vec3f> points(points_array[i].size());
    for(std::size_t j = 0; j < points.size(); ++j)
      points[j] = motions[i].second.transform(points_array[i][j]);

    saveOBJFile(out_filenameF2.c_str(), points, triangles_array[i]);
  }
  */
}

static void xml2tri(const std::string& in_filename, const std::string& out_filename_base)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;
  loadSceneFile(in_filename, points_array, triangles_array, motions);

  std::size_t n_obj = points_array.size();
  // save in local frame
  for(std::size_t i = 0; i < n_obj; ++i)
  {
    std::string out_filenameL = out_filename_base + num2string(i+1) + ".tri";
    savePolyDepthTriFile(out_filenameL.c_str(), points_array[i], triangles_array[i]);
  }
}


static void scenePenetrationTest(const std::string& filename, bool reverse = false, PenetrationDepthType pd_type = PDT_GENERAL_EULER)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;

  loadSceneFile(filename, points_array, triangles_array, motions);
  xml2tri(filename, "trimodel");

  if(reverse)
  {
    for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
    {
      Transform3f tf1 = motions[frame_id].first;
      Transform3f tf2 = motions[frame_id].second;
      motions[frame_id].first = tf2;
      motions[frame_id].second = tf1;
    }

    std::vector<Vec3f> p1 = points_array[0];
    std::vector<Vec3f> p2 = points_array[1];
    points_array[0] = p2;
    points_array[1] = p1;

    std::vector<Triangle> t1 = triangles_array[0];
    std::vector<Triangle> t2 = triangles_array[1];
    triangles_array[0] = t2;
    triangles_array[1] = t1;
  }


  std::ofstream motion_file("config.txt");
  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    Transform3f tf1 = motions[frame_id].first;
    Transform3f tf2 = motions[frame_id].second;

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf2.getRotation()(i, j) << " ";
    }
    motion_file << tf2.getTranslation()[0] << " " << tf2.getTranslation()[1] << " " << tf2.getTranslation()[2] << " ";
    motion_file << std::endl;

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf1.getRotation()(i, j) << " ";
    }
    motion_file << tf1.getTranslation()[0] << " " << tf1.getTranslation()[1] << " " << tf1.getTranslation()[2] << " ";
    motion_file << std::endl;


    motion_file << std::endl;

  }

  motion_file.close();

  

  BVHModel<OBBRSS>* m1 = new BVHModel<OBBRSS>();
  BVHModel<OBBRSS>* m2 = new BVHModel<OBBRSS>();
  m1->beginModel();
  m1->addSubModel(points_array[0], triangles_array[0]);
  m1->endModel();

  m2->beginModel();
  m2->addSubModel(points_array[1], triangles_array[1]);
  m2->endModel();

  Transform3f id;
  CollisionObject o1(boost::shared_ptr<CollisionGeometry>(m1), id);
  CollisionObject o2(boost::shared_ptr<CollisionGeometry>(m2), id);

  default_transform_distancer = DefaultTransformDistancer(o2.getCollisionGeometry());
  std::cout << "rotation weights ";
  default_transform_distancer.printRotWeight();
  std::size_t KNN_K = 10;
  LibSVMClassifier<6> classifier;
  
  std::vector<Transform3f> contact_vectors = penetrationDepthModelLearning(&o1, &o2, pd_type, &classifier, 1000000, 0, KNN_GNAT, KNN_K);

  classifier.save(filename + "model.txt");

  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    CollisionRequest c_request;
    CollisionResult c_result;
    collide(m1, motions[frame_id].first, m2, motions[frame_id].second, c_request, c_result);
    if(c_result.numContacts() == 0)
      std::cout << "0" << std::endl;
    else
    {
      PenetrationDepthRequest request(&classifier, default_transform_distance_func);
      request.contact_vectors = contact_vectors;
      PenetrationDepthResult result;

      penetrationDepth(m1, motions[frame_id].first, m2, motions[frame_id].second, request, result);
      
      std::cout << result.pd_value << " ";
      /*
      {
        Transform3f tf;
        relativeTransform2(result.resolved_tf, motions[frame_id].second, tf);
        std::cout << tf.getTranslation().length() << " ";
      }

      std::cout << (result.resolved_tf.getTranslation() - motions[frame_id].second.getTranslation()).length() << " ";
      */

      
      for(std::size_t i = 0; i < 3; ++i)
      {
        std::cout << result.resolved_tf.getTranslation()[i] << " ";
      }

      Matrix3f m = result.resolved_tf.getRotation();

      for(std::size_t i = 0; i < 3; ++i)
      {
        for(std::size_t j = 0; j < 3; ++j)
        {
          std::cout << m(i, j) << " ";
        }
      }
      std::cout << std::endl;
    }
    
  }
}



static void scenePenetrationTest_forBenchmark1a(const std::string& filename, PenetrationDepthType pd_type = PDT_GENERAL_EULER)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;

  loadSceneFile(filename, points_array, triangles_array, motions);
  xml2tri(filename, "trimodel");

  double delta_theta = boost::math::constants::pi<double>() / 10;


  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    Transform3f tf1 = motions[frame_id].first;
    Transform3f tf2 = motions[frame_id].second;
    double delta_z = tf1.getTranslation()[2] - motions[0].first.getTranslation()[2];

    motions[frame_id].first = motions[0].first;
    Matrix3f M;
    double theta = frame_id * delta_theta;
    M(0, 0) = cos(theta); M(0, 1) = -sin(theta); M(0, 2) = 0;
    M(1, 0) = sin(theta); M(1, 1) = cos(theta); M(1, 2) = 0;
    M(2, 0) = 0; M(2, 1) = 0; M(2, 2) = 1;
    tf2.setTransform(M, Vec3f(tf2.getTranslation()[0], tf2.getTranslation()[1], tf2.getTranslation()[2] - delta_z));
    motions[frame_id].second = tf2;
  }

  // now object 1 is fixed and object2 is moving

  std::ofstream motion_file("config.txt");
  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    Transform3f tf1 = motions[frame_id].first;
    Transform3f tf2 = motions[frame_id].second;

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf2.getRotation()(i, j) << " ";
    }
    motion_file << tf2.getTranslation()[0] << " " << tf2.getTranslation()[1] << " " << tf2.getTranslation()[2] << " ";
    motion_file << std::endl;

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf1.getRotation()(i, j) << " ";
    }
    motion_file << tf1.getTranslation()[0] << " " << tf1.getTranslation()[1] << " " << tf1.getTranslation()[2] << " ";
    motion_file << std::endl;


    motion_file << std::endl;

  }

  motion_file.close();

  

  BVHModel<OBBRSS>* m1 = new BVHModel<OBBRSS>();
  BVHModel<OBBRSS>* m2 = new BVHModel<OBBRSS>();
  m1->beginModel();
  m1->addSubModel(points_array[0], triangles_array[0]);
  m1->endModel();

  m2->beginModel();
  m2->addSubModel(points_array[1], triangles_array[1]);
  m2->endModel();

  Transform3f id;
  CollisionObject o1(boost::shared_ptr<CollisionGeometry>(m1), id);
  CollisionObject o2(boost::shared_ptr<CollisionGeometry>(m2), id);

  default_transform_distancer = DefaultTransformDistancer(o2.getCollisionGeometry());
  std::cout << "rotation weights ";
  default_transform_distancer.printRotWeight();
  std::size_t KNN_K = 10;
  LibSVMClassifier<6> classifier;
  
  std::vector<Transform3f> contact_vectors = penetrationDepthModelLearning(&o1, &o2, pd_type, &classifier, 10000, 0, KNN_GNAT, KNN_K);

  classifier.save(filename + "model.txt");

  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    CollisionRequest c_request;
    CollisionResult c_result;
    collide(m1, motions[frame_id].first, m2, motions[frame_id].second, c_request, c_result);
    if(c_result.numContacts() == 0)
      std::cout << "0" << std::endl;
    else
    {
      PenetrationDepthRequest request(&classifier, default_transform_distance_func);
      request.contact_vectors = contact_vectors;
      PenetrationDepthResult result;

      penetrationDepth(m1, motions[frame_id].first, m2, motions[frame_id].second, request, result);


      // std::cout << result.pd_value << " ";
      {
        Transform3f tf;
        relativeTransform2(result.resolved_tf, motions[frame_id].second, tf);
        std::cout << tf.getTranslation().length() << " ";
      }
      
      for(std::size_t i = 0; i < 3; ++i)
      {
        std::cout << result.resolved_tf.getTranslation()[i] << " ";
      }

      Matrix3f m = result.resolved_tf.getRotation();

      for(std::size_t i = 0; i < 3; ++i)
      {
        for(std::size_t j = 0; j < 3; ++j)
        {
          std::cout << m(i, j) << " ";
        }
      }
      std::cout << std::endl;
    }
    
  }
}



static void scenePenetrationR3Test(const std::string& filename, bool reverse = false)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;

  loadSceneFile(filename, points_array, triangles_array, motions);
  xml2tri(filename, "trimodel");

  // reverse when the first object is moving
  if(reverse)
  {
    for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
    {
      Transform3f tf1 = motions[frame_id].first;
      Transform3f tf2 = motions[frame_id].second;
      motions[frame_id].first = tf2;
      motions[frame_id].second = tf1;
    }

    std::vector<Vec3f> p1 = points_array[0];
    std::vector<Vec3f> p2 = points_array[1];
    points_array[0] = p2;
    points_array[1] = p1;

    std::vector<Triangle> t1 = triangles_array[0];
    std::vector<Triangle> t2 = triangles_array[1];
    triangles_array[0] = t2;
    triangles_array[1] = t1;
  }
  
  std::ofstream ground_truth_file("ground_truth.txt");
  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    Transform3f tf1 = motions[frame_id].first;
    Transform3f tf2 = motions[frame_id].second;

    /*
    double z1 = tf1.getTranslation()[2];
    double z2 = tf2.getTranslation()[2];

    double z = z1 - 0.25 - (1.5 + z2);
    if(z < 0)
      ground_truth_file << frame_id << " " << fabs(z) << std::endl;
    else
      ground_truth_file << frame_id << " " << 0 << std::endl;
    
    */
    
    double z1 = tf1.getTranslation()[2];
    double z2 = tf2.getTranslation()[2];

    double pd_z = fabs(z2 - 0.25 - (1.5 + z1));


    double x1 = tf1.getTranslation()[0];
    double x2 = tf2.getTranslation()[0];

    double pd_x = 0;
    if(x2 - 0.5 < x1 + 5)
    {
      pd_x = fabs(x2 - 0.5 - x1 - 5);
    }


    if(pd_x > 0)
    {
      ground_truth_file << frame_id << " " << std::min(pd_x, pd_z) << std::endl;
    }
    else
      ground_truth_file << frame_id << " " << 0 << std::endl;

    

    
  }
  ground_truth_file.close();


  return;

  std::ofstream motion_file("config.txt");
  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    Transform3f tf1 = motions[frame_id].first;
    Transform3f tf2 = motions[frame_id].second;

    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf2.getRotation()(i, j) << " ";
    }
    motion_file << tf2.getTranslation()[0] << " " << tf2.getTranslation()[1] << " " << tf2.getTranslation()[2] << " ";
    motion_file << std::endl;


    for(std::size_t i = 0; i < 3; ++i)
    {
      for(std::size_t j = 0; j < 3; ++j)
        motion_file << tf1.getRotation()(i, j) << " ";
    }
    motion_file << tf1.getTranslation()[0] << " " << tf1.getTranslation()[1] << " " << tf1.getTranslation()[2] << " ";
    motion_file << std::endl;

    motion_file << std::endl;
  }

  motion_file.close();

  

  BVHModel<OBBRSS>* m1 = new BVHModel<OBBRSS>();
  BVHModel<OBBRSS>* m2 = new BVHModel<OBBRSS>();
  m1->beginModel();
  m1->addSubModel(points_array[0], triangles_array[0]);
  m1->endModel();

  m2->beginModel();
  m2->addSubModel(points_array[1], triangles_array[1]);
  m2->endModel();

  Transform3f id;
  CollisionObject o1(boost::shared_ptr<CollisionGeometry>(m1), id);
  CollisionObject o2(boost::shared_ptr<CollisionGeometry>(m2), id);

  default_transform_distancer = DefaultTransformDistancer();
  default_transform_distancer.rot_x_weight = 0;
  default_transform_distancer.rot_y_weight = 0;
  default_transform_distancer.rot_z_weight = 0;
  std::cout << "rotation weights ";
  default_transform_distancer.printRotWeight();
  std::size_t KNN_K = 10;
  LibSVMClassifier<3> classifier;
  
  std::vector<Transform3f> contact_vectors = penetrationDepthModelLearning(&o1, &o2, PDT_TRANSLATIONAL, &classifier, 10000, 1, KNN_GNAT, KNN_K);

  classifier.save(filename + "model.txt");

  

  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    CollisionRequest c_request;
    CollisionResult c_result;
    collide(m1, motions[frame_id].first, m2, motions[frame_id].second, c_request, c_result);
    if(c_result.numContacts() == 0)
      std::cout << "0" << std::endl;
    else
    {
      PenetrationDepthRequest request(&classifier, default_transform_distance_func);
      request.contact_vectors = contact_vectors;
      PenetrationDepthResult result;

      penetrationDepth(m1, motions[frame_id].first, m2, motions[frame_id].second, request, result);


      std::cout << result.pd_value << " ";      
      
      for(std::size_t i = 0; i < 3; ++i)
      {
        std::cout << result.resolved_tf.getTranslation()[i] << " ";
      }

      Matrix3f m = result.resolved_tf.getRotation();

      for(std::size_t i = 0; i < 3; ++i)
      {
        for(std::size_t j = 0; j < 3; ++j)
        {
          std::cout << m(i, j) << " ";
        }
      }
      std::cout << std::endl;
    }
    
  }
}


static void sceneCollisionTest(const std::string& filename)
{
  std::vector<std::vector<Vec3f> > points_array;
  std::vector<std::vector<Triangle> > triangles_array;
  std::vector<std::pair<Transform3f, Transform3f> > motions;

  loadSceneFile(filename, points_array, triangles_array, motions);

  BVHModel<OBBRSS>* m1 = new BVHModel<OBBRSS>();
  BVHModel<OBBRSS>* m2 = new BVHModel<OBBRSS>();
  m1->beginModel();
  m1->addSubModel(points_array[0], triangles_array[0]);
  m1->endModel();

  m2->beginModel();
  m2->addSubModel(points_array[1], triangles_array[1]);
  m2->endModel();

  Transform3f id;
  CollisionObject o1(boost::shared_ptr<CollisionGeometry>(m1), id);
  CollisionObject o2(boost::shared_ptr<CollisionGeometry>(m2), id);

  for(std::size_t i = 0; i < motions.size(); ++i)
  {
    CollisionRequest request(100000, true);
    CollisionResult result;
    std::vector<Contact> contacts;
    collide(m1, motions[i].first, m2, motions[i].second, request, result);
    result.getContacts(contacts);
    std::cout << i << ": " << result.numContacts() << std::endl;
    if(contacts.size() > 0)
    {
      for(std::size_t j = 0; j < result.numContacts(); ++j)
        std::cout << "(" << contacts[j].b1 << "," << contacts[j].b2 << ")";
      std::cout << std::endl;
    }
  }
}


BOOST_AUTO_TEST_CASE(scene_test_penetration)
{
  boost::filesystem::path path(TEST_RESOURCES_DIR);

  
  RNG::setSeed(1);

  //std::cout << "manyframes/Model_7a" << std::endl;
  //std::string filename0 = (path / "manyframes/Model_7a.xml").string();
  //scenePenetrationTest_forBenchmark1a(filename0);

  
  //std::cout << "manyframes/Model_7a" << std::endl;
  //std::string filename0 = (path / "manyframes/Model_7a.xml").string();
  //scenePenetrationR3Test(filename0, true);
  //scenePenetrationTest(filename0, true);

  //std::cout << "manyframes/Model_7b" << std::endl;
  //std::string filename0 = (path / "manyframes/Model_7b.xml").string();
  //scenePenetrationR3Test(filename0, false);
  //scenePenetrationTest(filename0, false);

  
  std::cout << "manyframes/Model_6" << std::endl;
  std::string filename0 = (path / "manyframes/Model_6.xml").string();
  scenePenetrationTest(filename0, true);

  return;



  std::cout << "scenario-1-2-3/Model_1_Scenario_1" << std::endl;
  std::string filename1 = (path / "scenario-1-2-3/Model_1_Scenario_1.txt").string();
  scenePenetrationTest(filename1);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_1_Scenario_2" << std::endl;
  std::string filename2 = (path / "scenario-1-2-3/Model_1_Scenario_2.txt").string();
  scenePenetrationTest(filename2);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_1_Scenario_3" << std::endl;
  std::string filename3 = (path / "scenario-1-2-3/Model_1_Scenario_3.txt").string();
  scenePenetrationTest(filename3);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_2_Scenario_1" << std::endl;
  std::string filename4 = (path / "scenario-1-2-3/Model_2_Scenario_1.txt").string();
  scenePenetrationTest(filename4);

  std::cout << std::endl;
    
  std::cout << "scenario-1-2-3/Model_2_Scenario_2" << std::endl;
  std::string filename5 = (path / "scenario-1-2-3/Model_2_Scenario_2.txt").string();
  scenePenetrationTest(filename5);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_2_Scenario_3" << std::endl;
  std::string filename6 = (path / "scenario-1-2-3/Model_2_Scenario_3.txt").string();
  scenePenetrationTest(filename6);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_3_Scenario_1" << std::endl;
  std::string filename7 = (path / "scenario-1-2-3/Model_3_Scenario_1.txt").string();
  scenePenetrationTest(filename7);

  std::cout << std::endl;
    
  std::cout << "scenario-1-2-3/Model_3_Scenario_2" << std::endl;
  std::string filename8 = (path / "scenario-1-2-3/Model_3_Scenario_2.txt").string();
  scenePenetrationTest(filename8);

  std::cout << std::endl;

  std::cout << "scenario-1-2-3/Model_3_Scenario_3" << std::endl;
  std::string filename9 = (path / "scenario-1-2-3/Model_3_Scenario_3.txt").string();
  scenePenetrationTest(filename9);

  std::cout << std::endl;

  std::cout << "newdata/Model_1_Scenario_3" << std::endl;
  std::string filename10 = (path / "newdata/Model_1_Scenario_3.xml").string();
  scenePenetrationTest(filename10);

  std::cout << std::endl;


  std::cout << "newdata/Model_2_Scenario_3" << std::endl;
  std::string filename11 = (path / "newdata/Model_2_Scenario_3.xml").string();
  scenePenetrationTest(filename11);

  std::cout << std::endl;


  std::cout << "newdata/Model_3_Scenario_3" << std::endl;
  std::string filename12 = (path / "newdata/Model_3_Scenario_3.xml").string();
  scenePenetrationTest(filename12);

  std::cout << std::endl;

}

BOOST_AUTO_TEST_CASE(xml2obj_test)
{
  boost::filesystem::path path(TEST_RESOURCES_DIR);

  std::string filename_manyframe1 = (path / "manyframes/Model_1.xml").string();
  xml2obj(filename_manyframe1, "Model_1");

  return;

  std::string filename_manyframe0 = (path / "manyframes/Model_5.xml").string();
  xml2obj(filename_manyframe0, "Model_5");
  
  std::string filename_manyframe2 = (path / "manyframes/Model_4.xml").string();
  xml2obj(filename_manyframe2, "Model_4");
  
  std::string filename1 = (path / "scenario-1-2-3/Model_1_Scenario_1.txt").string();
  xml2obj(filename1, "Model_1_Scenario_1");

  std::string filename2 = (path / "scenario-1-2-3/Model_1_Scenario_2.txt").string();
  xml2obj(filename2, "Model_1_Scenario_2");
  
  std::string filename3 = (path / "scenario-1-2-3/Model_1_Scenario_3.txt").string();
  xml2obj(filename3, "Model_1_Scenario_3");

  std::string filename4 = (path / "scenario-1-2-3/Model_2_Scenario_1.txt").string();
  xml2obj(filename4, "Model_2_Scenario_1");

  std::string filename5 = (path / "scenario-1-2-3/Model_2_Scenario_2.txt").string();
  xml2obj(filename5, "Model_2_Scenario_2");

  std::string filename6 = (path / "scenario-1-2-3/Model_2_Scenario_3.txt").string();
  xml2obj(filename6, "Model_2_Scenario_3");

  std::string filename7 = (path / "scenario-1-2-3/Model_3_Scenario_1.txt").string();
  xml2obj(filename7, "Model_3_Scenario_1");

  std::string filename8 = (path / "scenario-1-2-3/Model_3_Scenario_2.txt").string();
  xml2obj(filename8, "Model_3_Scenario_2");

  std::string filename9 = (path / "scenario-1-2-3/Model_3_Scenario_3.txt").string();
  xml2obj(filename9, "Model_3_Scenario_3");

}
