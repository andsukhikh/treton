#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <cmath>

class Mesh
{
private:
	size_t PartitionNumber;
	std::vector<double> node_coord_;
public:
	Mesh(size_t size) : PartitionNumber(size)
	{
		node_coord_.reserve(PartitionNumber);
	}

	void init_mesh(double radius_first_boundary, double radius_second_boundary)
	{
		node_coord_[0]						= radius_first_boundary;
		node_coord_[PartitionNumber - 1]	= radius_second_boundary;

		//try to create own iterator which point to inner nodes
		double zero_layer = std::pow(radius_first_boundary, 2);
		double last_layer = std::pow(radius_second_boundary, 2);
		double elemntary_equal_area = (last_layer - zero_layer) / (2 * (PartitionNumber));

		for (size_t index = 1; index < PartitionNumber - 1; ++index)
		{
			double first_layer_radius = std::sqrt(zero_layer + (2 * index  + 1) * elemntary_equal_area);
			double second_layer_radius = std::sqrt(zero_layer + (2 * index + 3) * elemntary_equal_area);

			node_coord_[index] = (first_layer_radius + second_layer_radius) / 2;
		}
	}

	double& operator[] (int index)
	{
		return node_coord_[index];
	}
};

#endif 
