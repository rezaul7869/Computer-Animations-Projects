#define HEADERSIZE 54

#include "texture.h"


	
void swap(unsigned char* data, int l1, int l2)
{
	unsigned char temp = data[l1];
	data[l1] = data[l2];
	data[l2] = temp;
}

void convertFromBGRToRGB(unsigned char* data, int w, int h)
{
	for(int i = 0; i < (w*h*3); i+=3)
	{
		swap(data, i, i+2);
	}
}

GLuint loadBMP(const char* filename)
{
	unsigned char header[HEADERSIZE];
	unsigned int  pos;  
	unsigned int  w = 0;
	unsigned int  h = 0;
	unsigned int  size; 
	GLuint			textureId;
	unsigned char * data = NULL;
	if(data == NULL)
	{
		FILE* file = fopen(filename, "rb");

		if(!file)
			return -1;

		if(fread(header, 1, HEADERSIZE, file) != HEADERSIZE)
		{
			return -1;
		}

		if(header[0] != 'B' && header[1] != 'M')
		{
			return 0;
		}

		pos = *(int *) & header[0x0A];
		w = *(int *) & header[0x12];
		h = *(int *) & header[0x16];
		size = *(int *) & header[0x22];

		size = w * h * 3;
		if(pos == 0)	pos = 54;

		data = new unsigned char[size];
		if(fread(data, 1, size, file) != size)
			return -1;
		fclose(file);

		convertFromBGRToRGB(data, w, h);
	}

	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, w,h, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);


	if(data)
		delete(data);
	
	return textureId;
}