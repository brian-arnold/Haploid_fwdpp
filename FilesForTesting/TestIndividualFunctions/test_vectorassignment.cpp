#include <iostream>
#include <vector>

using namespace std ;

int main(){

	vector<int> vec ;
	vec.reserve(10) ;
	//vec[0] = 100 ;
	vec.push_back(100);
	for (int i=0; i<vec.size(); i++){
		cout << vec[i] << "\n" ;
	}
	
	return 0 ;
}
