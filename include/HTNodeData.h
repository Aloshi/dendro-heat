/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef INCLUDE_HTNODEDATA_H_
#define INCLUDE_HTNODEDATA_H_

#include <exception>

enum NodeDataIndices : int {
  U = 0,
  U_PRE = 1,
  HTNODEDATA_MAX = 2
};

/**
 * This class stores the data values for a single node.
 *
 * There are 3 items stored:
 * 1) the u value which is the calculated heat at the given point for the
 *    current time
 * 2) the u_pre value is the calculated heat for the previous time step
 * 3) u_analytical is the value of the analytical solution at the current time
 */
class HTNodeData {
 public:
  double u;
  double u_pre;

  /**
   * Returns reference to the given value in the object
   *
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double& value(int index) {
    switch (index) {
      case U: return u;
      case U_PRE: return u_pre;
      default: throw std::runtime_error("Invalid HTNodeData index");
    }
  }

  /**
   * Const reference version of value().
   * This function is required to be able to get read-only access to values
   * (e.g. when using a `const HTNodeData&` pointer or reference).
   * It is identical to the other value() function except for return type.
   *
   * @param index the index of the desired item
   * @returns const reference to the desired data item
   */
  const double& value(int index) const {
    switch (index) {
      case U: return u;
      case U_PRE: return u_pre;
      default: throw std::runtime_error("Invalid HTNodeData index");
    }
  }

  /**
   * Returns the name of the given data value in the object
   *
   * @param index the index of the desired item nsame
   * @return name of the specified data item
   */
  static char* name(int index) {
    static char str[256];
    if (index == U) {
      snprintf(str, sizeof(str), "u");
    } else if (index == U_PRE) {
      snprintf(str, sizeof(str), "u_pre");
    } else {
      throw std::runtime_error("Invalid HTNodeData index");
    }
    return str;
  }

  /**
   * Returns the number of the data items in the object
   *
   * @return number of the data items in the object
   */
  static int valueno() {
    return HTNODEDATA_MAX;
  }

  void UpdateDataStructures() {}
};

#endif  // INCLUDE_HTNODEDATA_H_
