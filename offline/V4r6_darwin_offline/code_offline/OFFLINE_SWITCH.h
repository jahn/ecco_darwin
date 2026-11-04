C     !ROUTINE: OFFLINE_SWITCH.h
C -------------------------------
C   OFFLINE_SWITCH.h
C  variable for switching on/off some calculations
C -------------------------------

C     offlineLoadGMRedi :: load from file GMRedi tensor (do not compute it)
C     offlineLoadKPP    :: load from file KPP mixing coeff (do not compute it)
C     offlineLoadConvec :: load from file Convective mixing (do not compute it)
C     offlineLoadAddMass:: load from file fluid source (do not compute it)
      COMMON /OFFLINE_SWITCH_L/
     &       offlineLoadGMRedi, offlineLoadKPP, offlineLoadConvec,
     &       offlineLoadAddMass
      LOGICAL offlineLoadGMRedi
      LOGICAL offlineLoadKPP
      LOGICAL offlineLoadConvec
      LOGICAL offlineLoadAddMass

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
