#ifndef __DataTypes_h__
#define __DataTypes_h__

typedef unsigned char	BOOLEAN;
//Char
typedef unsigned char	U8Bit; //(0 to 255)
typedef char			S8Bit; //(-128 to 127)

//short int
typedef unsigned short	U16Bit; //(0 to 65,535)
typedef signed short	S16Bit; //(-32,768 to 32,767)

//long int
typedef unsigned long	U32Bit; //(0 to 4,294,967,295)
typedef long			S32Bit; //(-2,147,483,648 to 2,147,483,647)

//float
typedef float			FLOAT4;

//double
typedef double			FLOAT8;

#define Z_TRUE			1
#define Z_FALSE			0

#define Z_SUCCESS		1
#define Z_FAILURE		0

#define Z_ON			1
#define Z_OFF			0

#define Z_ACTIVE		1
#define Z_INACTIVE		0

#define Z_CHANGE  		1
#define Z_NO_CHANGE  	0

#define Z_YES           1
#define Z_NO            0

#endif // __DataTypes_h__
