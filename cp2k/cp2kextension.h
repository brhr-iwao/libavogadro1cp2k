/**********************************************************************
  Cp2kExtension - CP2k Extension

  Copyright (C) Aoyama Iwao

  This file is not yet a part of the Avogadro molecular editor project.
  For more information about Avogadro, see <http://avogadro.cc/>

  Some code is based on Open Babel
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef CP2KEXTENSION_H
#define CP2KEXTENSION_H

#include <avogadro/extension.h>
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

namespace Avogadro {

  class Cp2kExtension : public Extension
  {
    Q_OBJECT
      AVOGADRO_EXTENSION("Cp2k", tr("Cp2k extension"),
                         tr("provides CHARMM style psf/pot files ganerator, Cp2k input generator and output reader "))

    public:
      // Constructor
      Cp2kExtension(QObject *parent=0);

      // Deconstructor
      virtual ~Cp2kExtension();

	  // This tells Avogadro what actions to create
      virtual QList<QAction *> actions() const;

	  // This returns a string that tells Avogadro where to put the menu entries
      virtual QString menuPath(QAction *action) const;

      virtual QDockWidget * dockWidget();

	  // This is called whenever a new molecule is loaded. 
      virtual QUndoCommand* performAction(QAction *action, GLWidget *widget);

     // This is called whenever a new molecule is loaded. 
      virtual void setMolecule(Molecule *molecule);

    private:
     // List of actions implemented by the extension
      QList<QAction *> m_actions;
      Molecule *m_molecule;

    private Q_SLOTS:

  };

  class Cp2kExtensionFactory : public QObject, public PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
    AVOGADRO_EXTENSION_FACTORY(Cp2kExtension)
  };

} // end namespace Avogadro

#endif
