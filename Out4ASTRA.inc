 990     format(1A6,2x,1(f8.3),2x,1(F8.3))
 991  format(1A6,2x,1I6)
        i=5 !output step
      nbnd2=0
      do j=1,NA1,i
      nbnd2=nbnd2+1
      enddo
      nbnd2=nbnd2+1
         open(20,file='IMAS.exp')
 write(20,*) ' iterIMAS'
  write(20,*) 'Name	Time	Value	Error'
         write(20,990) 'ABC', time, abc
         write(20,990) 'BTOR',time,BTOR
         write(20,990) 'RTOR',time,rtor
         write(20,990)  'AB', time, abc*1.1
         write(20,990) 'ELON',time, elon
         write(20,990) 'TRIAN',time,trian
         write(20,990) 'ELONG',time, elon
         write(20,990) 'TRICH',time,trian
         write(20,990) 'IPL',time,IPL
         write(20,990) 'AMJ',time,AMAIN(1)
         write(20,990) 'ZMJ',time,ZMAIN(1)
         write(20,991) 'NB1',nbnd2
         write(20,991) 'POINTS',nbnd2
          write(20,*) 'POINTS',nbnd2,'NAMEXP ','ne',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (ne(j),j=1,NA1,i),ne(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','ni',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (ni(j),j=1,NA1,i),ni(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','te',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (te(j),j=1,NA1,i),te(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','ti',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (ti(j),j=1,NA1,i),ti(NA1)
            write(20,*) 'POINTS',nbnd2,'NAMEXP ','zef',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (zef(j),j=1,NA1,i),zef(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','cu',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (cu(j),j=1,NA1,i),cu(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','mu',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (mu(j),j=1,NA1,i),mu(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','car5',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (fp(j),j=1,NA1,i),fp(NA1)
                   write(20,*) 'POINTS',nbnd2,'NAMEXP ','ametr',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (ametr(j),j=1,NA1,i),ametr(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','shif',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (shif(j),j=1,NA1,i),shif(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','g11',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (g11(j),j=1,NA1,i),g11(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','g22',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (g22(j),j=1,NA1,i),g22(NA1)
                   write(20,*) 'POINTS',nbnd2,'NAMEXP ','g33',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (g33(j),j=1,NA1,i),g33(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','ipol',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (ipol(j),j=1,NA1,i),ipol(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','pe',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (pe(j),j=1,NA1,i),pe(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','pi',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (pi(j),j=1,NA1,i),pi(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','car1',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (upl(j),j=1,NA1,i),upl(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','car4',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (cubs(j),j=1,NA1,i),cubs(NA1)
                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','car2',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (eqff(j),j=1,NA1,i),eqff(NA1)
                   write(20,*) 'POINTS',nbnd2,'NAMEXP ','car3',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (eqpf(j),j=1,NA1,i),eqpf(NA1)
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car6',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (cc(j),j=1,NA1,i),cc(NA1)
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car7',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J)),j=1,NA1,i),&
         FTLLMRy(BTOR,BMAXT(NA1),BDB0(NA1),BDB02(NA1),FOFB(NA1))
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car8',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (BMAXT(J),j=1,NA1,i), BMAXT(NA1)
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car9',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (BDB0(J),j=1,NA1,i),BDB0(NA1)
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car10',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (BDB02(J),j=1,NA1,i), BDB02(NA1)
                    write(20,*) 'POINTS',nbnd2,'NAMEXP ','car11',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
         write(20,*) time
         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
         write(20,*) (FOFB(J),j=1,NA1,i), FOFB(NA1)

!                  write(20,*) 'POINTS',nbnd2,'NAMEXP ','',' NTIMES 1 GRIDTY 12 FILTER 1.e-5 FACTOR 1.'
!         write(20,*) time
!         write(20,*) (XCP(j),j=1,NA1,i),XCP(NA1)
!         write(20,*) ((j),j=1,NA1,i),(NA1)
! FTLLMRy(BTOR,BMAXT(J),BDB0(J),BDB02(J),FOFB(J))
         close(20)
