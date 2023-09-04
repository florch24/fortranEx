      PROGRAM LANE_EMDEN

      implicit none  
      !--------------------------------------------------------------
       double precision m,a,x2,x1,h,b, fi
       complex*16 tita
       logical key
      !-------------------------------------------------------------- 
       open(10,file='condicion_inicial.dat',status='old') !no incluye número de ecuaciones
       read(10,*) fi, tita  

       open(11,file='indice_politropa.dat', status='old')
       read(11,*) m

       open(12,file='intervalo_integ.dat',status='old')
       read(12,*) a,b ! extremos de intervalo
    
       !paso de integración 
        h=0.1
       !corroboro lectura de datos
       !write(*,*) a, b, m, fi, tita

    !------------------Archivo de salida  --------------------------- 
        open(13,file='politropa.sal')
        write(13,*)'# epsilon  tita'
    !----------------------------------------------------------------
       x1=a
       key = .TRUE.
       DO WHILE(key)
           x2=DMIN1(x1+h,b)  
           CALL RK4(x2,tita,fi,h,m)
           write(13,*) x2,dble(tita)
           
           IF (.NOT.(b.GT.x2).OR. DABS(x2-b).LT. 0.5d-5)THEN
               key=.FALSE.
           END IF 
           x1=x2    

        END DO

      close(10)
      close(11)
      close(12)
      close(13) 

      END

!---------------------------------------------------------------------------------------
      SUBROUTINE RK4(x,y,z,h,m)
      !----------------------------------------------------  
      ! x corresponde a algún valor del intervalo
      !  y es tita_i compleja
      ! z es el fi_i
      !h es el paso
      !m es el índice politrópico  
      !----------------------------------------------------  
      implicit none
      double precision x,z, f, g,m,h
      complex*16 y
      !---- variables internas-------
      double precision k1, k2, k3, k4
      double precision l1, l2, l3, l4
      
      !----construimos las k y l ----
      
      k1=h*f(x,z)
      l1= h*g(x,y,m)
       
      k2=h*f(x+h/2., z+l1/2.)
      l2=h*g(x+h/2.,y+k1/2.,m)

      k3=h*f(x+h/2,z+l2/2.)
      l3=h*g(x+h/2.,y+k2/2.,m)

      k4= h*f(x+h/2,z+l3/2.)
      l4= h*g(x+h/2.,y+k3/2.,m)

      ! valores que devuelve la subrutina
      y= DCMPLX(y)+(1/6.)*(k1+k4+2*(k2+k3))
      z=z+(1/6.)*(l1+l4+2*(l2+l3))

      return
      END
       
      !---------------------------------------------------------------------------------  
      !--------------Expresiones para f y g --------------
      ! t=x, tita=y, fi=z

      DOUBLE PRECISION FUNCTION f(x,z)
      !corresponde a la función igualada a dtita/ depsilon
      double precision x,z

      f=z/(x**2)
      return
      END
      !------------------------------------------
      DOUBLE PRECISION FUNCTION g(x,y,m)
      !corresponde a la función igualada a dfi/depsilon
      !m es el índice politrópico
      double precision x,m
      complex*16 y

      g= (-1)*(dble(y**m))*(x**2)
      return
      END
      