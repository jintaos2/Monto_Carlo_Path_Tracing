#ifndef HDR_H
#define HDR_H
/***********************************************************************************
	Created:	17:9:2002
	FileName: 	hdrloader.h
	Author:		Igor Kravtchenko
	
	Info:		Load HDR image and convert to a set of float32 RGB triplet.
************************************************************************************/

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <glm.hpp>

typedef unsigned char RGBE[4];
#define _R 0
#define _G 1
#define _B 2
#define _E 3

#define MINELEN 8      // minimum scanline length for encoding
#define MAXELEN 0x7fff // maximum scanline length for encoding

class HDRLoaderResult
{
public:
    int width, height;
    // each pixel takes 3 float32, each component can be of any value...
    float *cols;
};

class HDRLoader
{
public:
    HDRLoaderResult color;
    bool success = false;

    HDRLoader(const char *fileName)
    {
        success = load(fileName);
    }
    glm::vec3 GetColor(int n)
    {
        return glm::vec3(color.cols[n * 3], color.cols[n * 3 + 1], color.cols[n * 3 + 2]);
    }
    inline bool load(const char *fileName)
    {
        int i;
        char str[200];
        FILE *file;

        file = fopen(fileName, "rb");
        if (!file)
            return false;

        fread(str, 10, 1, file);
        if (memcmp(str, "#?RADIANCE", 10))
        {
            fclose(file);
            return false;
        }

        fseek(file, 1, SEEK_CUR);

        char cmd[200];
        i = 0;
        char c = 0, oldc;
        while (true)
        {
            oldc = c;
            c = fgetc(file);
            if (c == 0xa && oldc == 0xa)
                break;
            cmd[i++] = c;
        }

        char reso[200];
        i = 0;
        while (true)
        {
            c = fgetc(file);
            reso[i++] = c;
            if (c == 0xa)
                break;
        }

        int w, h;
        if (!sscanf(reso, "-Y %ld +X %ld", &h, &w))
        {
            fclose(file);
            return false;
        }

        color.width = w;
        color.height = h;

        float *cols = new float[w * h * 3];
        color.cols = cols;

        RGBE *scanline = new RGBE[w];
        if (!scanline)
        {
            fclose(file);
            return false;
        }

        // convert image
        for (int y = h - 1; y >= 0; y--)
        {
            if (decrunch(scanline, w, file) == false)
                break;
            workOnRGBE(scanline, w, cols);
            cols += w * 3;
        }

        delete[] scanline;
        fclose(file);

        return true;
    }
    inline float convertComponent(int expo, int val)
    {
        float v = val / 256.0f;
        float d = (float)pow(2, expo);
        return v * d;
    }

    inline void workOnRGBE(RGBE *scan, int len, float *cols)
    {
        while (len-- > 0)
        {
            int expo = scan[0][_E] - 128;
            cols[0] = convertComponent(expo, scan[0][_R]);
            cols[1] = convertComponent(expo, scan[0][_G]);
            cols[2] = convertComponent(expo, scan[0][_B]);
            cols += 3;
            scan++;
        }
    }

    inline bool decrunch(RGBE *scanline, int len, FILE *file)
    {
        int i, j;

        if (len < MINELEN || len > MAXELEN)
            return oldDecrunch(scanline, len, file);

        i = fgetc(file);
        if (i != 2)
        {
            fseek(file, -1, SEEK_CUR);
            return oldDecrunch(scanline, len, file);
        }

        scanline[0][_G] = fgetc(file);
        scanline[0][_B] = fgetc(file);
        i = fgetc(file);

        if (scanline[0][_G] != 2 || scanline[0][_B] & 128)
        {
            scanline[0][_R] = 2;
            scanline[0][_E] = i;
            return oldDecrunch(scanline + 1, len - 1, file);
        }

        // read each component
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < len;)
            {
                unsigned char code = fgetc(file);
                if (code > 128)
                { // run
                    code &= 127;
                    unsigned char val = fgetc(file);
                    while (code--)
                        scanline[j++][i] = val;
                }
                else
                { // non-run
                    while (code--)
                        scanline[j++][i] = fgetc(file);
                }
            }
        }

        return feof(file) ? false : true;
    }

    inline bool oldDecrunch(RGBE *scanline, int len, FILE *file)
    {
        int i;
        int rshift = 0;

        while (len > 0)
        {
            scanline[0][_R] = fgetc(file);
            scanline[0][_G] = fgetc(file);
            scanline[0][_B] = fgetc(file);
            scanline[0][_E] = fgetc(file);
            if (feof(file))
                return false;

            if (scanline[0][_R] == 1 &&
                scanline[0][_G] == 1 &&
                scanline[0][_B] == 1)
            {
                for (i = scanline[0][_E] << rshift; i > 0; i--)
                {
                    memcpy(&scanline[0][0], &scanline[-1][0], 4);
                    scanline++;
                    len--;
                }
                rshift += 8;
            }
            else
            {
                scanline++;
                len--;
                rshift = 0;
            }
        }
        return true;
    }
};

#endif
