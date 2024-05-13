#include "MakeLibrary.h"
using namespace std;
namespace plt = matplot;

int main()
{
	plt::plot({ 1,2,3,4,5,6,7, }, "*");
	plt::show();
	cout << "Hello CMake." << endl;
	return 0;
}
