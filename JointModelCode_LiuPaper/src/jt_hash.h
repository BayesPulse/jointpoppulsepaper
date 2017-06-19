#pragma once
#include "jt_deconvolution_main.h"
Subject_type *initialize_subject(void);
Node_type *initialize_node(void);
Node_type *sequential_search(Node_type *list, double time);
void insert_node(Node_type *new_node, Node_type *list);
void insert_subject(Subject_type *new_subj, Subject_type *sublist);
void delete_node(Node_type *node, Node_type *list);
void print_list(Node_type *list);
void destroy_list(Node_type *list);
void destroy_sublist(Subject_type *sublist);