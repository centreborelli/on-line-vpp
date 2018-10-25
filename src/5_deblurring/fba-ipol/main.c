// Copyright (C) 2016, Jérémy Anger <jeremy.anger@cmla.ens-cachan.fr>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "vpp.h"
#include "fba.h"
#include "consistent_registration.h"

int main(int argc, char** argv)
{
    if (argc != 7) {
        fprintf(stderr, "usage:\n\t%s M W p downsampling input output \n", *argv);
        return EXIT_FAILURE;
    }

    int arg = 1;
    int M = atoi(argv[arg++]);
    int W = atoi(argv[arg++]);
    int p = atof(argv[arg++]);
    int downsampling = atoi(argv[arg++]);
    const char* input = argv[arg++];
    const char* output = argv[arg++];

    if (M < 1 || M > 10)
    {
        fprintf(stderr, "invalid 'M', should be a number between 1 and 10\n");
        return EXIT_FAILURE;
    }

    if (W < 2 || W % 2 != 0)
    {
        fprintf(stderr, "invalid 'W', should be an even number equal or greater than 2\n");
        return EXIT_FAILURE;
    }

    if (p > 40.f)
    {
        p = 40.f;
        fprintf(stderr, "'p' clamped to 40\n");
    }

    if (p < 0.f)
    {
        fprintf(stderr, "invalid 'p', should be a number between 0 and 40\n");
        return EXIT_FAILURE;
    }

    if (downsampling < 1 || downsampling > 10)
    {
        fprintf(stderr, "invalid 'downsampling', should be a number between 1 and 10\n");
        return EXIT_FAILURE;
    }

    int w, h, d;
    FILE* in = vpp_init_input(input, &w, &h, &d);
    assert(in);
    FILE* out = vpp_init_output(output, w, h, d);
    assert(out);

    image_float_t frame = new_image_float(w, h, d);
    image_float_t result = new_image_float(w, h, d);
    image_float_t buffer[M];
    image_float_t registered[M];
    for (int i = 0; i < M; i++) {
        buffer[i] = new_image_float(w, h, d);
        registered[i] = new_image_float(w, h, d);
    }

    int c = 0;
    int C = 0;
    while (vpp_read_frame(in, frame.data, w, h, d)) {
        if (C != M) C++;
        c = (c + 1) % C;
        copyf(buffer[c], frame);

        consistent_registration(registered, buffer, C, c, 1.f, 3);
        fba(result, registered, W, p, C);

        if (!vpp_write_frame(out, result.data, w, h, d))
            break;
    }

    return EXIT_SUCCESS;
}

