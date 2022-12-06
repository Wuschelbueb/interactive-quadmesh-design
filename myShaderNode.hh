//
// Created by wuschelbueb on 21.11.22.
//

#ifndef OPENFLIPPER_MYSHADERNODE_HH
#define OPENFLIPPER_MYSHADERNODE_HH

/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/






//=============================================================================
//
//  CLASS ShaderNode
//
//=============================================================================


//== INCLUDES =================================================================

#include <string>

#include <ACG/ShaderUtils/GLSLShader.hh>
#include <ACG/Scenegraph/BaseNode.hh>
#include <ACG/Scenegraph/DrawModes.hh>


/** \class ShaderNode ShaderNode.hh <ACG/Scenegraph/ShaderNode.hh>

    Set shaders for this node and all its children.
    All changes will be done in the enter() method undonecd ..cd
    in the leave() method.
**/

        class myShaderNode : public ACG::SceneGraph::BaseNode
        {
        public:

            /// Default constructor. Applies all properties.
            myShaderNode( BaseNode*           _parent = 0,
                        const std::string&  _name = "<ShaderNode>" );

            /// Destructor.
            virtual ~myShaderNode();


//            ACG_CLASSNAME(myShaderNode);

            /// set shader
            void enter(ACG::GLState& /*_state*/, const ACG::SceneGraph::DrawModes::DrawMode& _drawmode) override;
            /// disable shader
            void leave(ACG::GLState& /*_state*/, const ACG::SceneGraph::DrawModes::DrawMode& _drawmode) override;

            /// set shader
            void enterPick(ACG::GLState& /*_state*/, ACG::SceneGraph::PickTarget _target, const ACG::SceneGraph::DrawModes::DrawMode& _drawmode) override;
            /// disable shader
            void leavePick(ACG::GLState& /*_state*/, ACG::SceneGraph::PickTarget _target, const ACG::SceneGraph::DrawModes::DrawMode& _drawmode) override;


            /// Sets the shader dir.
            void setShaderDir( std::string _shaderDir);

            std::string shaderDir() { return shaderDir_; };

            std::string vertexShaderName(ACG::SceneGraph::DrawModes::DrawMode _drawmode, bool _pick = false);

            std::string fragmentShaderName(ACG::SceneGraph::DrawModes::DrawMode _drawmode, bool _pick = false);

            /** Removes the shader for the given draw mode
             * @param _drawmode Set the drawmode for which the shader should be deactivated
            */
            void disableShader (ACG::SceneGraph::DrawModes::DrawMode _drawmode);

            /** Sets a Shader for the given draw mode
             * @param _drawmode           Set the drawmode for which the shader should be activated
             * @param _vertexShader       filename of the Vertex Shader within the shader directory
             * @param _fragmentShader     filename of the Fragment Shader within the shader directory
             * @param _pickVertexShader   Vertex shader during picking
             * @param _pickFragmentShader Fragment shader during picking
            */

            void setShader( ACG::SceneGraph::DrawModes::DrawMode _drawmode ,
                            const std::string&  _vertexShader,
                            const std::string&  _fragmentShader,
                            std::string         _pickVertexShader = "",
                            std::string         _pickFragmentShader = "");

            /// Get the shader for the given drawMode
            GLSL::PtrProgram getShader( ACG::SceneGraph::DrawModes::DrawMode _drawmode, bool _pick = false);

            /// Check if a shader is available for the given drawMode
            bool hasShader( ACG::SceneGraph::DrawModes::DrawMode _drawmode, bool _pick = false);

            /// return available draw modes
            ACG::SceneGraph::DrawModes::DrawMode availableDrawModes() const override;

        private :
            // Path to the shaders ( if set ). If empty shaders will not be used.
            std::string shaderDir_;

            class ShaderInfo {

            public :
                ShaderInfo() :
                        vertexShader(0),
                        vertexShaderFile(""),
                        fragmentShader(0),
                        fragmentShaderFile(""),
                        program(0),
                        initialized(false)
                {
                };

                GLSL::PtrVertexShader   vertexShader;
                std::string             vertexShaderFile;

                GLSL::PtrFragmentShader fragmentShader;
                std::string             fragmentShaderFile;

                GLSL::PtrProgram        program;

                bool                    initialized;
            };

            std::map< size_t, ShaderInfo> shaders;
            std::map< size_t, ShaderInfo> pickShaders;
        };

#endif //OPENFLIPPER_MYSHADERNODE_HH
