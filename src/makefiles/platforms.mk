#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Felix
#
# Richard Beanland, Keith Evans & Rudolf A Roemer
#
# (C) 2013-17, all rights reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ifeq ($(MYPLATFORM),INT32GNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32GNU
endif

ifeq ($(MYPLATFORM),INT64NGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64NGNU
endif

ifeq ($(MYPLATFORM),INT64YGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64YGNU
endif

ifeq ($(MYPLATFORM),OPT32GNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32GNU
endif

ifeq ($(MYPLATFORM),OPT64NGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64NGNU
endif

ifeq ($(MYPLATFORM),OPT64YGNU)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64YGNU
endif

ifeq ($(MYPLATFORM),ITA64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeITA64Nifort
endif

ifeq ($(MYPLATFORM),ITA64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeITA64Yifort
endif

ifeq ($(MYPLATFORM),INT32ifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32ifort
endif

ifeq ($(MYPLATFORM),INT64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Nifort
endif

ifeq ($(MYPLATFORM),INT64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Yifort
endif

ifeq ($(MYPLATFORM),OPT32ifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32ifort
endif

ifeq ($(MYPLATFORM),OPT64Nifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Nifort
endif

ifeq ($(MYPLATFORM),OPT64Yifort)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Yifort
endif

ifeq ($(MYPLATFORM),INT32pgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT32pgf
endif

ifeq ($(MYPLATFORM),INT64Npgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Npgf
endif

ifeq ($(MYPLATFORM),INT64Ypgf)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Ypgf
endif

ifeq ($(MYPLATFORM),OPT32pgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT32pgf
endif

ifeq ($(MYPLATFORM),OPT64Npgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Npgf
endif

ifeq ($(MYPLATFORM),OPT64Ypgf)
include $(MYSTARTDIR)/makefiles/makefile.includeOPT64Ypgf
endif

ifeq ($(MYPLATFORM),AIXxlf90)
include $(MYSTARTDIR)/makefiles/makefile.includeAIXxlf90
endif

ifeq ($(MYPLATFORM),alpha)
include $(MYSTARTDIR)/makefiles/makefile.includealpha
endif

ifeq ($(MYPLATFORM),sun)
include $(MYSTARTDIR)/makefiles/makefile.includesun
endif

ifeq ($(MYPLATFORM),hpux)
include $(MYSTARTDIR)/makefiles/makefile.includehpux
endif

ifeq ($(MYPLATFORM),irix)
include $(MYSTARTDIR)/makefiles/makefile.includeirix
endif

ifeq ($(MYPLATFORM),sgi)
include $(MYSTARTDIR)/makefiles/makefile.includesgi
endif

ifeq ($(MYPLATFORM),INT64Nftn)
include $(MYSTARTDIR)/makefiles/makefile.includeINT64Nftn
endif
