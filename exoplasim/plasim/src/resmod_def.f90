      module resmod ! will be overwritten by MoSt

      parameter(NLAT_ATM = 32)
      parameter(NLEV_ATM = 10)
      parameter(NPRO_ATM = 1)
      parameter(NLIGHTS = 1)
      parameter(NSTEPS = 1)
      end module resmod

      
      ! T85L30 on 16/32/64 processors
      !parameter(NLAT_ATM = 128)
      !parameter(NLEV_ATM = 30)
      !!parameter(NPRO_ATM = 16)
      !parameter(NPRO_ATM = 32)
      !!parameter(NPRO_ATM = 64)


      ! T127L30 on 16/32/64 processors
      !parameter(NLAT_ATM = 192)
      !parameter(NLEV_ATM = 30)
      !!parameter(NPRO_ATM = 32)
      !parameter(NPRO_ATM = 48) ! Does not work so well, 32 more efficient? \




      ! T170L30 on 16/32/64 processors
      !parameter(NLAT_ATM = 256)
      !parameter(NLEV_ATM = 30)
      !parameter(NPRO_ATM = 32)
      !parameter(NPRO_ATM = 64) 