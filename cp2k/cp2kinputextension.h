/**********************************************************************
  Cp2kInputExtension - CP2k Input Extension

  Copyright (C) aoyama iwao

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.cc/>

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

#ifndef CP2KINPUTEXTENSION_H
#define CP2KINPUTEXTENSION_H

#include <avogadro/extension.h>
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

namespace Avogadro {

  class Cp2kInputExtension : public Extension
  {
    Q_OBJECT
      AVOGADRO_EXTENSION("Cp2k", tr("Cp2k input file extension"),
                         tr("provides cp2k input deck and psf/pot files for FIST calculation"))

    public:
      // Constructor
      Cp2kInputExtension(QObject *parent=0);

      // Deconstructor
      virtual ~Cp2kInputExtension();

	  // This tells Avogadro what actions to create
      virtual QList<QAction *> actions() const;

	  // This returns a string that tells Avogadro where to put the menu entries
	  // if branched menus (menuPath) are not necessary, this function will be unused.
      // virtual QString menuPath(QAction *action) const;

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

  class Cp2kInputExtensionFactory : public QObject, public PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
    AVOGADRO_EXTENSION_FACTORY(Cp2kInputExtension)
  };

} // end namespace Avogadro

#endif
