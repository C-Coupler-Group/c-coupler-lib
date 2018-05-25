/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PARSE_SPECIAL_WORDS
#define PARSE_SPECIAL_WORDS

#define SEPERATE_WORD_SPACE                      " "
#define SEPERATE_WORD_TAB                        "\t"
#define RESERVED_WORD_QUOTE_MARK                 "\""
#define RESERVED_WORD_EQUAL                      "="
#define RESERVED_WORD_LEFT_BRACKET               "("
#define RESERVED_WORD_RIGHT_BRACKET              ")"
#define RESERVED_WORD_COMMA                      ","
#define RESERVED_WORD_ATTRIBUTE                  "%"
#define RESERVED_WORD_ANNOTATION                 ';'

#define FUNCTION_WORD_ADD_NC_FILE                "add_nc_file"
#define FUNCTION_WORD_ADD_BIN_FILE               "add_bin_file"
#define FUNCTION_WORD_NEW_1D_GRID                "new_1D_grid"
#define FUNCTION_WORD_NEW_PARTIAL_GRID           "new_partial_grid"
#define FUNCTION_WORD_ADD_GRID_AREA              "add_grid_area"
#define FUNCTION_WORD_ADD_AREA_BOUND             "add_area_bound"
#define FUNCTION_WORD_COMBINE_GRIDS              "combine_grids"
#define FUNCTION_WORD_NEW_OPERATOR               "new_operator"
#define FUNCTION_WORD_COMBINE_OPERATORS          "combine_operators"
#define FUNCTION_WORD_WRITE_REMAP_WEIGHTS        "write_remap_wgts"
#define FUNCTION_WORD_READ_REMAP_WEIGHTS         "read_remap_wgts"
#define FUNCTION_WORD_REMAP                      "remap"
#define FUNCTION_WORD_WRITE_FIELD                "write_field"
#define FUNCTION_WORD_READ_FIELD                 "read_field"
#define FUNCTION_WORD_SET_BOUNDARY               "set_boundary"
#define FUNCIION_WORD_LEV_COORD_FROM_SIGMA       "lev_coord_from_sigma"
#define FUNCIION_WORD_LEV_COORD_FROM_HYBRID      "lev_coord_from_hybrid"
#define FUNCIION_WORD_SET_LEV_GRID_SIGMA_INFO    "set_lev_grid_sigma_info"
#define FUNCTION_WORD_READ_DATA                  "read_data"
#define FUNCTION_WORD_ALLOC_FIELD                "alloc_field"
#define FUNCTION_WORD_ISPAN                      "ispan"
#define FUNCTION_WORD_FSPAN                      "fspan"
#define FUNCTION_WORD_GEN_TEST_DATA              "generate_test_data"
#define FUNCTION_WORD_EVALUATE_ERROR             "evaluate_error"
#define FUNCTION_WORD_COMPUTE_REMAP_WEIGHTS      "calc_remap_wgts"
#define FUNCTION_WORD_EXTRACT_MASK               "extract_mask"
#define FUNCTION_WORD_COMPUTE_OCN_MASK           "calc_ocn_mask"
#define FUNCTION_WORD_SET_OPERATOR_PARA          "set_parameter"

#define OBJECT_TYPE_IO                           "IO_file"
#define OBJECT_TYPE_GRID                         "grid"
#define OBJECT_TYPE_FIELD_DATA                   "field_data"
#define OBJECT_TYPE_REMAP_OPERATOR               "remap_operator"
#define OBJECT_TYPE_REMAP_STRATEGY               "remap_scheme"
#define OBJECT_TYPE_REMAP_WEIGHTS                "remap_weights"

#endif
