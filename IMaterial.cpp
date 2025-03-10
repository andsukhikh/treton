#include "IMaterial.hpp"

IMaterial::IMaterial(MaterialParams&& params)
						: params_(std::move(params)) {}
