/** 
* \class Vector
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/



#ifndef VECTOR_H
#define VECTOR_H

/** \brief Vector_base provides the base for the class Vector by taking care of memory allocation
*
* The template struct provides the base for Vector class by allocating, reallocating 
* and freeing the memory dynamically in good old, fast and reliable C fashion.
* It is not recommendend to use the struct directly rather use the Vector class.
* 
* \param *_start_address ptr pointing to the address of the first element
* \param _size number of elements that can be stored
*
*/
template<typename T> struct Vector_base
{
 public:
        //std cto allocating space for exactly one element of type T
        inline Vector_base();
        //std cto allocating space for a number of elements of type T
        inline Vector_base(const int&);
        //std cto allocating space for a number of elements of type T and assigning a value direclty.
        inline Vector_base(const int&,const T&);
        //std dto
        inline ~Vector_base();
        //Allocates space for a number of elements 
        bool inline resize(const int&);
        //Allocates space for a number of elements and assigns a value to them
        bool inline resize(const int&,const T&);
                
        //points to the address of the first element
        T* _start_address;
        //number of elements that can be stored
        int _size;

};

/** \brief std cto allocating space for exactly one element of type T
*/
template<typename T> inline Vector_base<T>::Vector_base():_size(1)
{
 _start_address = (T*) malloc(sizeof (T));
}

/** \brief std cto allocating space for a number of elements of type T
*
* \param size is the number of elements space has to be allocated for
*/
template<typename T> inline Vector_base<T>::Vector_base(const int& size):_size(size)
{
 _start_address = (T*) malloc(sizeof (T) * size);
}

/** \brief std cto allocating space for a number of elements of type T and assigning 
* a value direclty.
*
* \param size is the number of elements space has to be allocated for
* \param value that shall be assigned to all elements
*/
template<typename T> inline Vector_base<T>::Vector_base(const int& size,const T& value):_size(size)
{  
  _start_address = NULL;  
  resize(size,value);  
}

/** \brief std dto
*/
template<typename T> inline Vector_base<T>::~Vector_base()
{
}

/** \brief Allocates space for a number of elements
*
* \param new_size the new amount of elements that can be stored
*/
template<typename T> bool inline Vector_base<T>::resize(const int& new_size)
{
 T* tmp = (T*)realloc(_start_address,sizeof(T) * new_size);
 if (tmp != NULL) {_start_address = tmp;_size = new_size;return true; }
 else { return false; }
}

/** \brief Allocates space for a number of elements and assigns a value to them
*
* \param new_size the new amount of elements that can be stored
* \param value is assigned to the elements
*/
template<typename T> bool inline Vector_base<T>::resize(const int& new_size,const T& value)
{
 T* tmp = (T*) realloc(_start_address,sizeof(T) * new_size);
 if (tmp != NULL) 
 {
   _start_address = tmp;
   _size = new_size;
    
       
   if (this->_size > 1 )
   for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
   {
      *((this->_start_address)+j)=value;
      *((this->_start_address)+j+1)=value;
   }
   if (this->_size%2==1)
       *((this->_start_address)+this->_size-1)=value;

   return true;
   
   
 }   
 else{ return false; }
}

/** \brief Specialiced resize for float that allocates space for a number of float elements 
* and assigns a value to them
*
* \param new_size the new amount of elements that can be stored
* \param value is assigned to the elements
*/
template<> bool inline Vector_base<float>::resize(const int& new_size,const float& value)
{
 float* tmp = (float*) realloc(_start_address,sizeof(float) * new_size);
 if (tmp != NULL) 
 {
    _start_address = tmp;
    _size = new_size;
    
   int             typesize = sizeof(float);
   const float*    leftv    = _start_address;     
   int             vecsize  = this->_size;
   const float     scalar   = (float)value;
      
   asm(  
       "push esi;"
       "push edi;"
       "push ebx;"
       "mov edi, %0;"
       "mov esi, %1;"
       "mov ebx, %2;"
       "mov ecx, %3;"
       "loopb:;"
       "fld DWORD PTR [esi];"      
       "fstp DWORD PTR [edi];"       
       "add edi, ebx;"
       "dec ecx;"
       "jz toendb;"
       "fld DWORD PTR [esi];"      
       "fstp DWORD PTR [edi];"       
       "add edi, ebx;"
       "dec ecx;"
       "jnz loopb;"
       "toendb:;"
       "pop ebx;"
       "pop edi;"
       "pop esi;" : : "r"(leftv),"r"(&scalar),"r"(typesize) , "r"(vecsize) 
       );        
   return true;

   
 }   
 else{ return false; }
}

/**
* @brief The template class implements a mathemamatical n-dim vector of any type, which 
* has also the ability to be used as a container.
*
* The Vector class is build upon the struct Vector_base implementing a n-dim
* mathematical vector of any type extending the struct by coming up with the standard
* operations of adding and subtracting vectors as well as multiplying a vector by a scalar.
* The class permits the multiplication of a vector of arbitrary type e.g. int by a scalar
* of an arbitrary type e.g. float in order to allow arbitrary multiplications.
* Thus arbitrary multiplications are possible where an internal cast to thr type of the 
* vector is done.
* The class is optimized by assembler progamming for the datatype FLOAT. It is suggested to use
* instead of int also the FLOAT datatype, since it does run at least two times faster than any 
* optimized pure C++ code.
* Since the class can be used as container the operators +,- and = applied to the container are
* passed to the operators of the contained correspondng data type.
*
* Usage: 
*        Vector<float> v1(2,0.1); defining a 2-dim vector with the value 0.1 in each component
*        Vector<float> v2(2,0.4); defining a 2-dim vector with the value 0.4 in each component
*        int scalar = 2;          
*        v1+= v2;                 adding vector v1 and v2 where the result is assigned to the non-const v1
*        v1 = v1 * scalar;        multiplying vector v1 by the scalar and assigning it to vector v1
*
* \param _dummy dummy variable for operator[] that is returned if given index is out of range
*/



template<typename T> class Vector : protected Vector_base<T>
{
 typedef Vector_base<T>    _base;
 
 public: 
         // std cto
         inline Vector();
         // cto defining a vector with given dimension 
         inline Vector(const int&);
         // cto defining a vector with given dimension and values
         inline Vector(const int&,const T&);
         // std dto
         inline ~Vector();
         //operator= assigns a vector to the lhs
         inline Vector<T>& operator=(const Vector<T>&);
         //operator+= based on the operations = and +
         inline Vector<T>& operator+=(const Vector<T>&);
         //operator-= based on the operations = and -
         inline Vector<T>& operator-=(const Vector<T>&);
         //operator*= for scalars based on the operations = and *
         template<typename S> inline Vector<T>& operator*=(const S&);
         //operator/= for scalars based on the operations = and /
         template<typename S> inline Vector<T>& operator/=(const S&);
         // operator+ adding two vectors and returning a new result vector         
         inline Vector<T> operator+(const Vector<T>&);
         // const operator+ adding two const vectors and returning a new result vector         
         inline Vector<T> operator+(const Vector<T>&) const;
         // operator- subtracting two const vectors and returning a new result vector
         inline Vector<T> operator-(const Vector<T>&);
         // operator- subtracting two const vectors and returning a new result vector
         inline Vector<T> operator-(const Vector<T>&) const;
         // operator* multiplies a vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<T> operator*(const S&);
         // operator* const multiplies a vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<T> operator*(const S&) const;
         // operator/ divides a vector by an arbitrary scalar from the right        
         template<typename S> inline Vector<T> operator/(const S&);
         // operator/ const divides a vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<T> operator/(const S&) const;
         // operator[] returns a reference to the indexed element  
         inline T& operator[](const int&);
         // operator[] const returns a reference to the indexed element  
         inline const T& operator[](const int&) const;
         /** \brief friend operator* permitting a multiplication from the left with an arbitrary type
         *
         * Operators where the left element is not a class element have to be global or friend operators.
         * This friend operator* permits a multiplication from the left with an arbitray type 
         * where the scalar is internaly casted to the type of the vector.
         * \param factor is the scalar by which the vector is multiplied
         * \param v is the vector that is multiplied by the scalar 
         */   
         template<typename S> friend inline Vector<T> operator*(const S& factor,const Vector<T>& v)
         {
          return v*factor;
         }
         // changes the vector's dimension
         inline bool resize(const int&);
         // changes the vector큦 dimension and assigns a value to each entry
         inline bool resize(const int&, const T&);
         // size returns the dimension (number of components) of the vector 
         inline int size() const;
 private:
         // dummy variable for operator[] that is returned if given index is out of range
         T _dummy;         
};

/**
* @brief The specialized class implements the mathemamatical n-dim FLOAT vector 
* with very fast operations implemented in assembler.
*
* The Vector class is build upon the struct Vector_base implementing an optimized n-dim
* mathematical FLOAT vector extending the struct by coming up with the standard
* operations of adding and subtracting vectors as well as multiplying a vector by a scalar
* of any time.
* The class permits the multiplication of a FLOAT vector by a scalar an arbitrary type 
* e.g. int , double in order to allow arbitrary multiplications.
* Thus arbitrary multiplications are possible where an internal cast to thr type of the 
* vector is done.
* The class is optimized by assembler progamming. It is suggested to use the class as well
* for int, since it does run at least two times faster than any optimized pure C++ code.
*
* Usage: 
*        Vector<float> v1(2,0.1); defining a 2-dim vector with the value 0.1 in each component.
*        Vector<float> v2(2,0.4); defining a 2-dim vector with the value 0.4 in each component.
*        int scalar = 2;          
*        v1+= v2;                 adding vector v1 and v2 where the result is assigned to the non-const v1.
*        v1 = v1 * scalar;        multiplying vector v1 by the scalar and assigning it to vector v1.
*
* \param _dummy dummy variable for operator[] that is returned if given index is out of range
*/

template<> class Vector<float> : protected Vector_base<float>
{
 typedef Vector_base<float>    _base;
 
 public: 
         //std cto
         inline Vector();
         //cto defining a vector with given dimension 
         inline Vector(const int&);
         //cto defining a vector with given dimension and values
         inline Vector(const int&,const float&);
         //std destructor
         inline ~Vector();
         //operator= assigns a vector to the lhs
         inline Vector<float>& operator=(const Vector<float>&);
         //operator+= based on the operations = and +
         inline Vector<float>& operator+=(const Vector<float>&);
         //operator-= based on the operations = and -
         inline Vector<float>& operator-=(const Vector<float>&);
         //operator*= for scalars based on the operations = and *
         template<typename S> inline Vector<float>& operator*=(const S&);
         //operator/= for scalars based on the operations = and /
         template<typename S> inline Vector<float>& operator/=(const S&);
         //operator+ adding two const FLOAT vectors and returning a new FLOAT result vector         
         inline Vector<float> operator+(const Vector<float>&);
         //operator+ adding two FLOAT vectors and returning a new FLOAT result vector
         inline Vector<float> operator+(const Vector<float>&) const;
         //operator- subtracting two const FLOAT vectors and returning a new FLOAT result vector
         inline Vector<float> operator-(const Vector<float>&);
         //const operator- subtracting two const FLOAT vectors and returning a new FLOAT result vector
         inline Vector<float> operator-(const Vector<float>&) const;
         //operator* multiplies a FLOAT vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<float> operator*(const S&);
         //operator* const multiplies a FLOAT vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<float> operator*(const S&) const;
         //operator/ divides a FLOAT vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<float> operator/(const S&);
         //operator/ const divides a vector by an arbitrary scalar from the right 
         template<typename S> inline Vector<float> operator/(const S&) const;
         // operator[] returns a reference to the indexed element  
         inline float& operator[](const int&);
         // operator[] const returns a reference to the indexed element  
         inline const float& operator[](const int&) const;
         /** \brief friend operator* permitting a multiplication from the left with an arbitrary type
         *
         * Operators where the left element is not a class element have to be global or friend operators.
         * This friend operator* permits a multiplication from the left with an arbitray type 
         * where the scalar is internaly casted to the type of the vector.
         * \param factor is the scalar by which the vector is multiplied
         * \param v is the vector that is multiplied by the scalar 
         */
         template<typename S> friend inline Vector<float> operator*(const S& factor,const Vector<float>& v)
         {
          return v*factor;
         }
         //changes the vector's dimension
         inline bool resize(const int&);
         //changes the vector큦 dimension and assigns a value to each entry
         inline bool resize(const int&, const float&);
         //returns the dimension (number of components) of the vector 
         inline int size() const;
 private:
         //dummy variable for operator[] that is returned if given index is out of range
         float _dummy;
             
};


/** @brief std cto
*/
template<typename T> inline Vector<T>::Vector():_base(){}

/** @brief std cto
*/
inline Vector<float>::Vector():_base(){}

/** @brief cto defining a vector with given dimension 
*   @param size gives the number of coordinates equivalent to the dimension of the vector
*/
template<typename T> inline Vector<T>::Vector(const int& size):_base(size){}

/** @brief cto defining a vector with given dimension 
*   @param size gives the number of coordinates equivalent to the dimension of the vector
*/
inline Vector<float>::Vector(const int& size):_base(size){}

/** @brief cto defining a vector with given dimension and values
*   @param size gives the number of coordinate entries equivalent to the dimension of the vector
*   @param value that is assignend to each dimension entry
*/
template<typename T> inline Vector<T>::Vector(const int& size,const T& value):_base(size,value){}

/** @brief cto defining a vector with given dimension and values
*   @param size gives the number of coordinate entries equivalent to the dimension of the vector
*   @param value that is assignend to each dimension entry
*/
inline Vector<float>::Vector(const int& size,const float& value):_base(size,value){}

/** @brief std destructor
*/ 
template<typename T> inline Vector<T>::~Vector()
{free(this->_start_address);this->_start_address=NULL;}

/** @brief std destructor
*/ 
inline Vector<float>::~Vector(){
free(this->_start_address);
this->_start_address=NULL;
}


/** \brief operator= assigns a vector to the lhs
*
* The operator= assigns the rhs to the lhs where if the size of both does not match
* the lhs of the assignment is resized to the size of the rhs.
* \param to_assign is the rhs to be assigned to the lhs
*/
template<typename T> inline Vector<T>& Vector<T>::operator=(const Vector<T>& to_assign)
{
  if ( this != &to_assign )

  {
    if ( this->_size != to_assign._size )
    {
      
      this->resize(to_assign._size);     
      
    }  
    if (this->_size > 1 )
      for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
      {
            *((this->_start_address)+j)= *(to_assign._start_address+j);
            *((this->_start_address)+j+1)= *(to_assign._start_address+j+1);
      }
    if (this->_size%2==1)
             *((this->_start_address)+this->_size-1)= *(to_assign._start_address+this->_size-1);       
 
  }
  return *this;  
}

/** \brief operator= assigns a vector to the lhs
*
* The operator= assigns the rhs to the lhs where if the size of both does not match
* the lhs of the assignment is resized to the size of the rhs.
* \param to_assign is the rhs to be assigned to the lhs
*/
inline Vector<float>& Vector<float>::operator=(const Vector<float>& to_assign)
{
  if ( this != &to_assign )

  {
    if ( this->_size != to_assign._size )
    {
      this->resize(to_assign._size);
    }  
    
    void*       leftv       =     &this->operator[](0);    
    const void* rightv      =     &to_assign[0];
    int         typesize    =     sizeof(float);
    int         vecsize     =     this->_size;
    
    asm(
       "push esi;"
       "push edi;"
       "push ebx;"
       "mov edi, %0;"
       "mov esi, %1;"
       "mov ebx, %2;"
       "mov ecx, %3;"
       "loopa:;"
       "fld DWORD PTR [esi];"
       "fstp DWORD PTR [edi];"
       "add esi, ebx;"
       "add edi, ebx;"
       "dec ecx;"
       "jz toenda;"
       "fld DWORD PTR [esi];"  
       "fstp DWORD PTR [edi];"
       "add esi, ebx;"
       "add edi, ebx;"
       "dec ecx;"
       "jnz loopa;"
       "toenda:;"
       "pop ebx;"
       "pop edi;"
       "pop esi;" : : "r"(leftv),"r"(rightv),"r"(typesize) , "r"(vecsize) 
       );
    
    this->_size          = to_assign._size;
   
  }  
  return *this; 
}  

/** operator+= based on the operations = and +
*
* \param to_add_assign is what has to be added
*/
template<typename T> inline Vector<T>& Vector<T>::operator+=(const Vector<T>& to_add_assign)
{
  (*this) = (*this) + to_add_assign;
  return *this;
}  

/** operator+= based on the operations = and +
*
* \param to_add_assign is what has to be added
*/
inline Vector<float>& Vector<float>::operator+=(const Vector<float>& to_add_assign)
{
  (*this) = (*this) + to_add_assign;
  return *this;
}  

/** operator-= based on the operations = and -
*
* \param to_subtract_assign is what has to be subtracted 
*/
template<typename T> inline Vector<T>& Vector<T>::operator-=(const Vector<T>& to_subtract_assign)
{
  (*this) = (*this) - to_subtract_assign;
  return *this;
}  

/** operator-= based on the operations = and -
*
* \param to_subtract_assign is what has to be subtracted 
*/
inline Vector<float>& Vector<float>::operator-=(const Vector<float>& to_subtract_assign)
{
  (*this) = (*this) - to_subtract_assign;
  return *this;
}  

/** operator*= for scalars based on the operations = and *
*
* \param to_mul_assign is the factor by which the vector has to be multiplied
*/
template<typename T> template<typename S> inline Vector<T>& Vector<T>::operator*=(const S& to_mul_assign)
{
  (*this) = (*this) * to_mul_assign;
  return *this;
}  

/** operator*= for scalars based on the operations = and *
*
* \param to_mul_assign is the factor by which the vector has to be multiplied
*/
template<typename S> inline Vector<float>& Vector<float>::operator*=(const S& to_mul_assign)
{
  (*this) = (*this) * to_mul_assign;
  return *this;
}  

/** operator/= for scalars based on the operations = and /
*
* \param to_div_assign is the factor by which the vector has to be divided
*/
template<typename T> template<typename S> inline Vector<T>& Vector<T>::operator/=(const S& to_div_assign)
{
  (*this) = (*this) / to_div_assign;
  return *this;
}  

/** operator/= for scalars based on the operations = and /
*
* \param to_div_assign is the factor by which the vector has to be divided
*/
template<typename S> inline Vector<float>& Vector<float>::operator/=(const S& to_div_assign)
{
  (*this) = (*this) / to_div_assign;
  return *this;
}  

/** \brief const operator+ adding two const vectors and returning a new result vector
*
*   The const operator+ adds two const vectors of an arbitray type, creates a result vector
*   and returns that sum vector.
*   \param to_add is the vector that is added to lhs element
*/
template<typename T> inline Vector<T> Vector<T>::operator+(const Vector<T>& to_add) const
{
   
 if ( this->_size == to_add._size )
 {
     Vector<T>* result   =new Vector<T>(this->_size);
     if (this->_size > 1 )
          for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
          {
             (*result)[j]   =  *(to_add._start_address+j) + *((this->_start_address)+j) ;
             (*result)[j+1] = *((this->_start_address)+j+1) + *(to_add._start_address+j+1);
          }
       if (this->_size%2==1)
           (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) + *(to_add._start_address+this->_size-1);       
              
       return (*result);
 }  
} 

/** \brief const operator+ adding two const FLOAT vectors and returning a new FLOAT result vector
*
*   The const operator+ adds two const FLOAT vectors, creates a FLOAT result vector
*   and returns that sum vector.
*   \param to_add is the vector that is added to lhs element
*/
inline Vector<float> Vector<float>::operator+(const Vector<float>& to_add) const
{
 if ( this->_size == to_add._size )
 {
     Vector<float>* result   =new Vector<float>(this->_size);
 

     const float* leftv       =     &this->operator[](0);    
     float*       resultptr   =     &(*result)[0];
     int          typesize    =     sizeof(float);
     int          vecsize     =     this->_size;
     const float* rightv      =     &to_add[0];
     
     asm(  
         "push esi;"
         "push edi;"
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "looppc:;"
         "fld DWORD PTR [edi];"
         "fadd DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jz toendpc;"
         "fld DWORD PTR [edi];"
         "fadd DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jnz looppc;"
         "toendpc:;" 
         "pop ebx;"
         "pop edi;"
         "pop esi;": : "r"(leftv),"r"(rightv),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );
     return (*result);
  
 }
 return *this;
}

/** \brief operator+ adding two vectors and returning a new result vector
*
*   The const operator+ adds two vectors of an arbitray type, creates a result vector
*   and returns that sum vector.
*   \param to_add is the vector that is added to lhs element
*/
template<typename T> inline Vector<T> Vector<T>::operator+(const Vector<T>& to_add)
{

 Vector<T>* result   =  new Vector<T>(this->_size);
 
 if ( this->_size == to_add._size )
 {     
     if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   =  *(to_add._start_address+j) + *((this->_start_address)+j) ;
           (*result)[j+1] = *((this->_start_address)+j+1) + *(to_add._start_address+j+1);
        }
     if (this->_size%2==1)
         (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) + *(to_add._start_address+this->_size-1);       
            
     
  
 }
 return (*result);
}

/** \brief operator+ adding two FLOAT vectors and returning a new FLOAT result vector
*
*   The operator+ adds two FLOAT vectors, creates a FLOAT result vector
*   and returns that sum vector.
*   \param to_add is the vector that is added to lhs element
*/
inline Vector<float> Vector<float>::operator+(const Vector<float>& to_add)
{
  
 Vector<float>* result   =  new Vector<float>(this->_size); 
 
 if ( this->_size == to_add._size )
 {
     const float* leftv      = &this->operator[](0);    
     float*       resultptr  = &(*result)[0];
     int          typesize   = sizeof(float);
     int          vecsize    = this->_size;
     const float* rightv     = &to_add[0];
     
     asm(  
         "push esi;"
         "push edi;"
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loopp:;"
         "fld DWORD PTR [edi];"
         "fadd DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jz toendp;"
         "fld DWORD PTR [edi];"
         "fadd DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jnz loopp;"
         "toendp:;" 
         "pop ebx;"
         "pop edi;"
         "pop esi;": : "r"(leftv),"r"(rightv),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );     
 }
 return (*result); 

}

/** \brief const operator- subtracting two const vectors and returning a new result vector
*
*   The const operator- subtracts two const vectors of an arbitray type, creates a result vector
*   and returns that difference vector.
*   \param to_subtract is the vector that is subtracted from the lhs element
*/
template<typename T> inline Vector<T> Vector<T>::operator-(const Vector<T>& to_subtract) const
{ 
 if ( this->_size == to_subtract._size )
 { 
     Vector<T>* result   = new Vector<T>(this->_size);
     
     if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) - *(to_subtract._start_address+j)  ;
           (*result)[j+1] = *((this->_start_address)+j+1) - *(to_subtract._start_address+j+1);
        }
     if (this->_size%2==1)
         (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) - *(to_subtract._start_address+this->_size-1);       
            
     return (*result);
   
 }
}  

/** \brief const operator- subtracting two const FLOAT vectors and returning a new FLOAT result vector
*
*   The const operator- subtracts two const FLOAT vectors, creates a FLOAT result vector
*   and returns that difference vector.
*   \param to_subtract is the vector that is subtracted from the lhs element
*/
inline Vector<float> Vector<float>::operator-(const Vector<float>& to_subtract) const
{
  
 if ( this->_size == to_subtract._size )
 {
       
     Vector<float>* result    =   new Vector<float>(this->_size);
     const float*   leftv     =   &this->operator[](0);    
     float*         resultptr =   &(*result)[0];
     int            typesize  =   sizeof(float);
     int            vecsize   =   this->_size;
     const float*   rightv    =   &to_subtract[0];
     
     asm(  
         "push esi;"
         "push edi;"
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loopsc:;"
         "fld DWORD PTR [edi];"
         "fsub DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jz toendsc;"
         "fld DWORD PTR [edi];"
         "fsub DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jnz loopsc;"
         "toendsc:;"
         "pop ebx;"
         "pop edi;"
         "pop esi;" : : "r"(leftv),"r"(rightv),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );
     return (*result);
 }
 return *this;
}

/** \brief operator- subtracting two vectors and returning a new result vector
*
*   The operator- subtracts two vectors of an arbitray type, creates a result vector
*   and returns that difference vector.
*   \param to_subtract is the vector that is subtracted from the lhs element
*/
template<typename T> inline Vector<T> Vector<T>::operator-(const Vector<T>& to_subtract)
{ 
 
 Vector<T>* result   =  new Vector<T>(this->_size); 
 
 if ( this->_size == to_subtract._size )
 {
 
     if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) - *(to_subtract._start_address+j)  ;
           (*result)[j+1] = *((this->_start_address)+j+1) - *(to_subtract._start_address+j+1);
        }
     if (this->_size%2==1)
         (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) - *(to_subtract._start_address+this->_size-1);       
 
 }
 return (*result);
}  

/** \brief operator- subtracting two FLOAT vectors and returning a new FLOAT result vector
*
*   The operator- subtracts two FLOAT vectors, creates a FLOAT result vector
*   and returns that difference vector.
*   \param to_subtract is the vector that is subtracted from the lhs element
*/
inline Vector<float> Vector<float>::operator-(const Vector<float>& to_subtract)
{
 
 Vector<float>* result   =  new Vector<float>(this->_size); 
 
 if ( this->_size == to_subtract._size )
 {       
     const float*  leftv       =   &this->operator[](0);    
     float*        resultptr   =   &(*result)[0];
     int           typesize    =   sizeof(float);
     int           vecsize     =   this->_size;
     const float*  rightv      =   &to_subtract[0];
     
     asm(  
         "push esi;"
         "push edi;"
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loops:;"
         "fld DWORD PTR [edi];"
         "fsub DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jz toends;"
         "fld DWORD PTR [edi];"
         "fsub DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "add esi, ebx;"
         "dec ecx;"
         "jnz loops;"
         "toends:;"
         "pop ebx;"
         "pop edi;"
         "pop esi;" : : "r"(leftv),"r"(rightv),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );
 }
 return (*result);
}

/** \brief operator* multiplies a vector by an arbitrary scalar from the right 
*
* The operator* multiplies a vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to the type of the vector.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is multiplied
*/
template<typename T> template<typename S> inline Vector<T> Vector<T>::operator*(const S& factor)
{   

   Vector<T>* result   =   new Vector<T>(this->_size);
   const T    scalar   =   (T)factor;
   
   if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) * scalar   ;
           (*result)[j+1] = *((this->_start_address)+j+1) * scalar;
        }
   if (this->_size%2==1)
       (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) * scalar;       
   return (*result);

}

/** \brief operator* multiplies a FLOAT vector by an arbitrary scalar from the right 
*
* The operator* multiplies a FLOAT vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to float.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is multiplied
*/
template<typename S> inline Vector<float> Vector<float>::operator*(const S& factor)
{   

   Vector<float>*   result      =   new Vector<float>(this->_size);
   const float*     leftv       =   &this->operator[](0);    
   float*           resultptr   =   &(*result)[0];
   int              typesize    =   sizeof(float);
   int              vecsize     =   this->_size;
   const float      scalar      =   (float)factor;

   asm(  
       "push esi;"
       "push edi;"
       "push ebx;"
       "mov edi, %0;"
       "mov esi, %1;"
       "mov ebx, %2;"
       "mov ecx, %3;"
       "mov edx, %4;"
       "loopm:;"
       "fld DWORD PTR [edi];"
       "fmul DWORD PTR [esi];"
       "fstp DWORD PTR [edx];"       
       "add edi, ebx;"
       "add edx, ebx;"
       "dec ecx;"
       "jz toendm;"
       "fld DWORD PTR [edi];"
       "fmul DWORD PTR [esi];"
       "fstp DWORD PTR [edx];"       
       "add edi, ebx;"
       "add edx, ebx;"
       "dec ecx;"
       "jnz loopm;"
       "toendm:;"
       "pop ebx;"
       "pop edi;"
       "pop esi;" : : "r"(leftv),"r"(&scalar),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
          
         );

  return (*result);
}  

/** \brief operator* const multiplies a vector by an arbitrary scalar from the right 
*
* The operator* multiplies a vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to the type of the vector.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is multiplied
*/
template<typename T> template<typename S> inline Vector<T> Vector<T>::operator*(const S& factor) const
{   

   Vector<T>* result   =new Vector<T>(this->_size);
   const T scalar= (T)factor;
   if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) * scalar   ;
           (*result)[j+1] = *((this->_start_address)+j+1) * scalar;
        }
   if (this->_size%2==1)
       (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) * scalar;       
   return (*result);

}

/** \brief operator* const multiplies a FLOAT vector by an arbitrary scalar from the right 
*
* The operator* multiplies a FLOAT vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to float.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is multiplied
*/
template<typename S> inline Vector<float> Vector<float>::operator*(const S& factor) const
{
 
     Vector<float>* result      =    new Vector<float>(this->_size);
     const float*   leftv       =    &this->operator[](0);    
     float*         resultptr   =    &(*result)[0];
     int            typesize    =    sizeof(float);
     int            vecsize     =    this->_size;
     const float    scalar      =    (float)factor;

     asm(  
         "push esi;"
         "push edi;"
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loopmulc:;"
         "fld DWORD PTR [edi];"
         "fmul DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "dec ecx;"
         "jnz loopmulc;"
         "pop ebx;"
         "pop edi;"
         "pop esi;" : : "r"(leftv),"r"(&scalar),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );

  return (*result);
}  

/** \brief operator/ divides a vector by an arbitrary scalar from the right 
*
* The operator/ divides a vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to the type of the vector.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is divided
*/
template<typename T> template<typename S> inline Vector<T> Vector<T>::operator/(const S& factor) 
{   

   Vector<T>* result   =  new Vector<T>(this->_size);
   const T scalar      = (T)factor;
   if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) / scalar   ;
           (*result)[j+1] = *((this->_start_address)+j+1) / scalar;
        }
   if (this->_size%2==1)
       (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) / scalar;       
   return (*result);

}

/** \brief operator/ divides a FLOAT vector by an arbitrary scalar from the right 
*
* The operator* divides FLOAT vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to float.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is divided
*/
template<typename S> inline Vector<float> Vector<float>::operator/(const S& factor) 
{
  
     Vector<float>* result    =   new Vector<float>(this->_size);
     const float*   leftv     =   &this->operator[](0);    
     float*         resultptr =   &(*result)[0];
     int            typesize  =   sizeof(float);
     int            vecsize   =   this->_size;
     const float    scalar    =   (float)factor;
     
     asm(
         "push esi;"
         "push edi;"  
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loopd:;"
         "fld DWORD PTR [edi];"
         "fdiv DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "dec ecx;"
         "jz toendd;"
         "fld DWORD PTR [edi];"
         "fdiv DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "dec ecx;"
         "jnz loopd;"
         "toendd:;"
         "pop ebx;"
         "pop edi;"
         "pop esi;" : : "r"(leftv),"r"(&scalar),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );
  
  return (*result);
}  

/** \brief operator/ const divides a vector by an arbitrary scalar from the right 
*
* The operator* divides a vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator/.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to the type of the vector.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is divided
*/
template<typename T> template<typename S> inline Vector<T> Vector<T>::operator/(const S& factor) const
{   

   Vector<T>* result   =  new Vector<T>(this->_size);
   const T    scalar   =  (T)factor;
   if (this->_size > 1 )
        for (int i=0, j=0; i < (this->_size / 2) ; i++,j+=2)
        {
           (*result)[j]   = *((this->_start_address)+j) / scalar   ;
           (*result)[j+1] = *((this->_start_address)+j+1) / scalar;
        }
   if (this->_size%2==1)
       (*result)[this->_size-1] = *((this->_start_address)+this->_size-1) / scalar;       
   return (*result);

}

/** \brief operator/ const divides a FLOAT vector by an arbitrary scalar from the right 
*
* The operator* divides FLOAT vector by an arbitrary scalar from the right 
* and is therefore the complement of the friend operator*.
* If the scalar is of a different than the vector then the type of the scalar is internaly
* casted to float.
* The result is the lhs vector.
* \param factor is the scalar by which the vector is divided
*/
template<typename S> inline Vector<float> Vector<float>::operator/(const S& factor) const
{
  
     Vector<float>* result    =    new Vector<float>(this->_size);
     const float*   leftv     =    &this->operator[](0);    
     float*         resultptr =    &(*result)[0];
     int            typesize  =    sizeof(float);
     int            vecsize   =    this->_size;
     const float    scalar    =    (float)factor;
     
     asm(
         "push esi;"
         "push edi;"  
         "push ebx;"
         "mov edi, %0;"
         "mov esi, %1;"
         "mov ebx, %2;"
         "mov ecx, %3;"
         "mov edx, %4;"
         "loopdc:;"
         "fld DWORD PTR [edi];"
         "fdiv DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "dec ecx;"
         "jz toenddc;"
         "fld DWORD PTR [edi];"
         "fdiv DWORD PTR [esi];"
         "fstp DWORD PTR [edx];"       
         "add edi, ebx;"
         "add edx, ebx;"
         "dec ecx;"
         "jnz loopdc;"
         "toenddc:;"
         "pop ebx;"
         "pop edi;"
         "pop esi;" : : "r"(leftv),"r"(&scalar),"r"(typesize) , "r"(vecsize) , "r"(resultptr)
            
           );
  
  return (*result);
}  
/** \brief operator[] returns a reference to the indexed element  
* \param index is the coordinate of vector that is going to be accessed
*/
template<typename T> inline T& Vector<T>::operator[](const int& index)
{
 if ( 0 <= index && index < this->_size) 
    return *(this->_start_address + index);
  else
    return _dummy;
}

/** \brief operator[] returns a reference to the indexed element  
* \param index is the coordinate of vector that is going to be accessed
*/
inline float& Vector<float>::operator[](const int& index)
{
 if ( 0 <= index && index < this->_size) 
    return *(this->_start_address + index);
  else
    return _dummy;
}

/** \brief operator[] const returns a reference to the indexed element  
* \param index is the coordinate of vector that is going to be accessed
*/
template<typename T> const inline T& Vector<T>::operator[](const int& index) const
{
 if ( 0 <= index && index < this->_size) 
    return *(this->_start_address + index);
  else
    return _dummy;
}

/** \brief operator[] const returns a reference to the indexed element  
* \param index is the coordinate of vector that is going to be accessed
*/
const inline float& Vector<float>::operator[](const int& index) const
{
 if ( 0 <= index && index < this->_size) 
    return *(this->_start_address + index);
  else
    return _dummy;
}

/** \brief changes the vector's dimension
*   \param new_size is the new size of the vector
*/
template<typename T> inline bool Vector<T>::resize(const int& new_size)
{
 return (_base::resize(new_size));
}

/** \brief changes the vector's dimension
*   \param new_size is the new size of the vector
*/
inline bool Vector<float>::resize(const int& new_size)
{
 return (_base::resize(new_size));
}

/** \brief changes the vector큦 dimension and assigns a value to each entry
*   \param new_size is the new size of the vector
*   \param value is the value that is assigned to each vector entry
*/
template<typename T> inline bool Vector<T>::resize(const int& new_size, const T& value)
{
 return (_base::resize(new_size,value));
}

/** \brief changes the vector큦 dimension and assigns a value to each entry
*   \param new_size is the new size of the vector
*   \param value is the value that is assigned to each vector entry
*/
inline bool Vector<float>::resize(const int& new_size, const float& value)
{
 return (_base::resize(new_size,value));
}

/** \brief size returns the dimension (number of components) of the vector 
*/
template<typename T> inline int Vector<T>::size() const
{
 return this->_size;
}

/** \brief size returns the dimension (number of components) of the vector 
*/
inline int Vector<float>::size() const
{
 return this->_size;
}


#endif
