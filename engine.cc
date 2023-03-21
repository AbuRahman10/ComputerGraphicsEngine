#include "easy_image.h"
#include "ini_configuration.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include "Line2D.h"
#include "Functies.h"
#include "Figure.h"

#include "l_parser.h"
#include <fstream>

using namespace std;
using namespace img;
using namespace LParser;


EasyImage colorrectangle(unsigned int width, unsigned int height)
{
    img::EasyImage image(width, height);
    for (unsigned int i = 0; i < height; i++)
    {
        for (unsigned int j = 0; j < width; j++)
        {
            image(i, j).red = i;
            image(i, j).green = j;
            image(i, j).blue = (i + j) % 256;
        }
    }
    return image;
}

EasyImage blocks(unsigned int Wi, unsigned int Hi, unsigned int Nx, unsigned int Ny, vector<double> wit, vector<double> zwart, bool invertColors)
{
    int Wb = Wi / Nx;
    int Hb = Hi / Ny;
    img::EasyImage image(Wi, Hi);

    for (int x = 0; x < Hi; x++)
    {
        for (int y = 0; y < Wi; y++)
        {
            int Bx = x/Wb;
            int By = y/Hb;

            if (invertColors == false)
            {
                if ((Bx+By)%2 == 0)
                {
                    image(x,y).red = wit[0] * 255;
                    image(x,y).green = wit[1] * 255;
                    image(x,y).blue = wit[2] * 255;
                }
                else
                {
                    image(x,y).red = zwart[0] * 255;
                    image(x,y).green = zwart[1] * 255;
                    image(x,y).blue = zwart[2] * 255;
                }
            }
            else
            {
                if ((Bx+By)%2 == 0)
                {
                    image(x,y).red = zwart[0] * 255;
                    image(x,y).green = zwart[1] * 255;
                    image(x,y).blue = zwart[2] * 255;
                }
                else
                {
                    image(x,y).red = wit[0] * 255;
                    image(x,y).green = wit[1] * 255;
                    image(x,y).blue = wit[2] * 255;
                }
            }
        }
    }

    return image;
}

EasyImage quarterCircle(unsigned int Wi, unsigned int Hi, int N, vector<double> backgroundColor, vector<double> lineColor)
{
    int Hs = Hi / (N - 1); // verticale afstand
    int Ws = Wi / (N - 1); // horizontale afstand

    img::EasyImage image(Wi, Hi);

    img::Color color(lineColor[0] * 255, lineColor[1] * 255, lineColor[2] * 255);

    image.draw_line(0,Hi-1,Wi-1,Hi-1,color);
    for (int i = 0; i < Hi-1; i++)
    {
        if (Ws*i < Wi and Hs*i < Hi)
        {
            image.draw_line(0,Ws*i,Hs*i,Hi-1,color);
        }
    }

    return image;
}

EasyImage eye(unsigned int Wi, unsigned int Hi, int N, vector<double> backgroundColor, vector<double> lineColor)
{
    int Hs = Hi / (N - 1); // verticale afstand
    int Ws = Wi / (N - 1); // horizontale afstand

    img::EasyImage image(Wi, Hi);

    img::Color color(lineColor[0] * 255, lineColor[1] * 255, lineColor[2] * 255);

    image.draw_line(0,Hi-1,Wi-1,Hi-1,color);
    for (int i = 0; i < Hi-1; i++)
    {
        if (Ws*i < Wi and Hs*i < Hi)
        {
            image.draw_line(0,Ws*i,Hs*i,Hi-1,color);
        }
    }
    image.draw_line(0,0,Hi-1,0,color);

    for (int i = 0; i < Hi-1; i++)
    {
        if (Ws*i < Wi and Hs*i < Hi)
        {
            image.draw_line(Ws*i,0,Hi-1,Hs*i,color);
        }
    }
    image.draw_line(Hi-1,0,Hi-1,Hi-1,color);


    return image;
}

EasyImage diamond(unsigned int Wi, unsigned int Hi, int N, vector<double> backgroundColor, vector<double> lineColor)
{
    int Hs = Hi / (N - 1); // verticale afstand
    int Ws = Wi / (N - 1); // horizontale afstand

    img::EasyImage image(Wi, Hi);

    img::Color color(lineColor[0] * 255, lineColor[1] * 255, lineColor[2] * 255);

    image.draw_line(0,Hi-1,Wi-1,Hi-1,color);
    for (int i = 0; i < Hi-1; i++)
    {
        if (Ws*i < Wi and Hs*i < Hi)
        {
            image.draw_line(0,Ws*i,Hs*i,Hi-1,color);
        }
    }
    image.draw_line(0,0,Hi-1,0,color);

    for (int i = 0; i < Hi-1; i++)
    {
        if (Ws*i < Wi and Hs*i < Hi)
        {
            image.draw_line(Ws*i,0,Hi-1,Hs*i,color);
        }
    }
    image.draw_line(Hi-1,0,Hi-1,Hi-1,color);

    return image;
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    string type = configuration["General"]["type"].as_string_or_die();
    //string input_file = configuration["2DLSystem"]["inputfile"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> colors = configuration["Figure0"]["color"].as_double_tuple_or_die();
    vector<double> backgroundColors = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();

    Colour colour(colors[0],colors[1],colors[2]);

    if (type == "2DLSystem")
    {
        LSystem2D l_system;

        ifstream input_stream("input_file");
        input_stream >> l_system;
        input_stream.close();

        Lines2D lines2D = Functies::drawLSystem(l_system);
        return Functies::draw2DLines(lines2D, size, backgroundColors);
    }
    else if (type == "Wireframe")
    {
        Figures3D figures3D;
        // INLEZEN VAN DE INI FILE
        Lines2D lines2D = Functies::pasFigure(figures3D,configuration,colour);
        // TEKENEN VAN DE 2D LIJNEN
        return Functies::draw2DLines(lines2D, size, backgroundColors);
    }
    else
    {
        EasyImage image(500,500);
        return image;
    }
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
