///////////////////////////////////////////////
// README

// Copyright (c) 2024 Nils Jonas Norberg

// Inspired by the following document:
// (1980) Euclidean Distance Mapping, PER-ERIK DANIELSSON, Department of Electrical Engineering, Linkoping University

///////////////////////////////////////////////
// INTERFACE

#ifndef OK_SDF_H
#define OK_SDF_H

#include <cstdint> // for uint8_t

// takes a multi-channel image and changes it to a 1-channel mask with 0/1 per pixel
void ok_sdf_extract_mask_in_place(uint8_t* data_inout, int w, int h, int channels, int extract_channel, int threshold, bool invert);

// copies the src to the destination (dest has to be pre-allocated to fit the padded image)
void ok_sdf_pad_mask(uint8_t* dst_img, const uint8_t* src_img, int w, int h, int pad);

// takes a 1-channel mask and replaces it with a sdf-image
void ok_sdf_process_in_place(uint8_t* data_inout, int w, int h, int scale, int spread, int& w_out, int& h_out, uint8_t* rgba_to_ps_fix = nullptr, bool verbose = true);

// takes a 1-channel image and shrinks to fit (in-place)
void ok_sdf_crop_image_in_place(uint8_t* src_img, int& w_inout, int &h_inout);

#endif // OK_SDF_H
