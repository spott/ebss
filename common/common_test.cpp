#include<common/common.hpp>
#include<vector>
#include<iostream>

int
main (int argc, char** argv)
{
    std::vector<int> a;
    a.push_back(1);
    a.push_back(3);
    a.push_back(2);

    std::cout << "export the vector: " << std::endl;
    std::cout << a[0] << ", " << a[1] << ", " << a[2] << std::endl;
    common::export_vector_binary("fname", a);

    std::cout << "import the vector: " << std::endl;
    std::vector<int> c = common::import_vector_binary<int>("fname");
    std::cout << c[0] << ", " << c[1] << ", " << c[2] << std::endl;

    std::cout << "change the vector: " << std::endl;
    std::vector<float> b = common::vector_type_change<int,float>(c);  //This is what I want to work..
    std::cout << b[0] + .01 << ", " << b[1] + .01 << ", " << b[2] + .01 << std::endl;

    std::cout << "import then change vector: " << std::endl;
    std::vector<float> d = common::vector_type_change<int,float>(
            common::import_vector_binary<int>("fname")
            );  //This is what I want to work..
    std::cout << d[0] + .01 << ", " << d[1] + .01 << ", " << d[2] + .01 << std::endl;
}
