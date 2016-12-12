#include<iostream>
#include<fstream>
#include"tests.h"
//#include"gtest/gtest.h"


using namespace std;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}