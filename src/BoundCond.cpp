#include <vector>
#include <string>
#include <iostream>
extern int N, fict;
extern std::string x_left_bound, x_right_bound;

void BoundCond(std::vector<std::vector<double>>& W) {

	// Free bound
	if (x_left_bound == "free") {
		for (int i = 0; i < fict; i++){
			// Density
			/*W[fict - 1 - i][0] = W[fict + i][0];
			// Velosity
			W[fict - 1 - i][1] = W[fict + i][1];		
			// Pressure
			W[fict - 1 - i][2] = W[fict + i][2];*/
			
			W[fict - 1 - i][0] = W[fict][0];
			// Velosity
			W[fict - 1 - i][1] = W[fict][1];		
			// Pressure
			W[fict - 1 - i][2] = W[fict][2];
		}
	}
	if (x_right_bound == "free") {
		for (int i = 0; i < fict; i++){
			// Density
			/*W[N + fict - 1 + i][0] = W[N + fict - 2 - i][0];
			// Velosity
			W[N + 2*fict - 2 - i][1] = W[N + fict - 2][1];
			// Pressure
			W[N + 2*fict - 2 - i][2] = W[N + fict - 2][2];*/
			
			// Density
			W[N + fict - 1 + i][0] = W[N + fict - 2][0];
			// Velosity
			W[N + 2*fict - 2 - i][1] = W[N + fict - 2][1];
			// Pressure
			W[N + 2*fict - 2 - i][2] = W[N + fict - 2][2];
		}
	}

	/*// Wall
	if (bound_name == "wall") {
		for (int i = 0; i < fict; i++) {
			// Density
			W[i][0] = W[2 * fict - 1 - i][0];
			W[N - 1 - i][0] = W[N - 2 * fict + i][0];

			// Velosity
			W[i][1] = -W[2 * fict - 1 - i][1];
			W[N - 1 - i][1] = -W[N - 2 * fict - i][1];
			
			// Pressure
			W[i][2] = W[2 * fict - 1 - i][2];
			W[N - 1 - i][2] = W[N - 2 * fict + i][2];
		}
	}

	// Periodic
	if (bound_name == "period") {
		for (int i = 0; i < fict; i++) {
			// Deinsity
			W[i][0] = W[N - 2 * fict + i][0];
			W[N - fict + i][0] = W[fict + i][0];

			// Velosity
			W[i][1] = W[N - 2 * fict + i][1];
			W[N - fict + i][1] = W[fict + i][1];

			// Pressure
			W[i][2] = W[N - 2 * fict + i][2];
			W[N - fict + i][2] = W[fict + i][2];
		}
	}
	*/
	return;
}
