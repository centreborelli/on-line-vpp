/*
 *
 * This file is part of the Image Processing Framework.
 *
 * Copyright(c) 2011 Miguel Colom.
 * miguel.colom@cmla.ens-cachan.fr
 *
 * This file may be licensed under the terms of of the
 * GNU General Public License Version 2 (the ``GPL'').
 *
 * Software distributed under the License is distributed
 * on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied. See the GPL for the specific language
 * governing rights and limitations.
 *
 * You should have received a copy of the GPL along with this
 * program. If not, go to http://www.gnu.org/licenses/gpl.html
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include "CFramework.h"
#include "CInspector.h"
#include "../algo.h"

#include <time.h>

int main(int argc, char **argv) {
  CFramework *fw = CFramework::get_framework();  
  try {
    algorithm(argc, argv);
    delete fw;
  }
  catch(...) {    
    printf("A exception raised while running the algorithm.\n");
    printf("Backtrace:\n");
    fw->print_backtrace(stdout);
    // [ToDo]: call user handler for the exception
    exit(-1);
  }
}
