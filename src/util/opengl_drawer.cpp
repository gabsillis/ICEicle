#include <iceicle/opengl_drawer.hpp>


namespace iceicle::gl { 

    const char *arrow_vert_shader =
#include "../../../shaders/arrow2d_shader.vert"
    ;

    const char *arrow_frag_shader =
#include "../../../shaders/arrow2d_shader.frag"
    ;

    const char *arrow_geom_shader =
#include "../../../shaders/arrow2d_shader.geom"
    ;

    template<>
    ShapeDrawer<Arrow>::ShapeDrawer() : vao(), shader(arrow_vert_shader, arrow_frag_shader, arrow_geom_shader)
    {
		vertex_attributes.push_back(0);
		vertex_attributes.push_back(1);
		vao.bind();
		vao.buffers.emplace("arrow_data", GL_ARRAY_BUFFER);
    }

    template<>
    void ShapeDrawer<Arrow>::buffer_data(){
		vao.bind();
		vao["arrow_data"].bind();
		vao["arrow_data"].buffer_data(draw_list.size(), draw_list.data(), GL_STATIC_DRAW);

		vao.bind();
		vao["arrow_data"].set_attr_pointer<GLfloat>(0, 3, 6, 0);
		vao["arrow_data"].set_attr_pointer<GLfloat>(1, 3, 6, 3);
    }

    template<>
    void ShapeDrawer<Arrow>::draw_arrays(){
	    glDrawArrays(GL_POINTS, 0, 2);
    }

    // ===================
    // = Generated Arrow =
    // ===================

    const char *arrow2_vert = 
#include "../../shaders/bounding_box_scale.vert"
    ;

    const char *arrow2_frag = R"(
    #version 330 core 
    out vec3 color;
    void main(){
	    color = vec3(0.6, 0.0, 0.8);
    }
    )";

    template<>
    ShapeDrawer<ArrowGenerated>::ShapeDrawer() 
    : vao(), shader(arrow2_vert, arrow2_frag, nullptr)
    {
	    vertex_attributes.push_back(0);
	    vao.bind();
	    vao.buffers.emplace("tri_data", GL_ARRAY_BUFFER);
    }

    template<>
    void ShapeDrawer<ArrowGenerated>::buffer_data(){
	    vao.bind();
	    vao["tri_data"].bind();
	    vao["tri_data"].buffer_data<glm::vec3>(9 * draw_list.size(), draw_list[0].pts, GL_STATIC_DRAW);
	    vao["tri_data"].set_attr_pointer<GLfloat>(0, 3, 3, 0);
    }

    template<>
    void ShapeDrawer<ArrowGenerated>::draw_arrays(){
	    glDrawArrays(GL_TRIANGLES, 0, 9 * draw_list.size());
    }


    const char *tri_vert_shader = 
#include "../../shaders/bounding_box_scale.vert"
    ;
    const char *tri_frag_shader = 
#include "../../shaders/triangle2d.frag"
    ;
    template<>
    ShapeDrawer<Triangle>::ShapeDrawer() : vao(), shader(tri_vert_shader, tri_frag_shader, nullptr)
    {
	    vertex_attributes.push_back(0);
	    vao.bind();
	    vao.buffers.emplace("vertex_data", GL_ARRAY_BUFFER);
    }

    template<>
    void ShapeDrawer<Triangle>::buffer_data(){
	    vao.bind();
	    vao["vertex_data"].bind();
	    vao["vertex_data"].buffer_data<glm::vec3>(draw_list.size() * 3, &(draw_list.data()->pt1), GL_STATIC_DRAW);
	    vao["vertex_data"].set_attr_pointer<GLfloat>(0, 3, 3, 0);
    }

    template<>
    void ShapeDrawer<Triangle>::draw_arrays(){
	    glDrawArrays(GL_TRIANGLES, 0, 3 * draw_list.size());
    }


    // ===================
    // = Buffered Shapes =
    // ===================

    const char *curve_vert = 
#include "../../shaders/bounding_box_scale.vert"
    ;
    const char *curve_frag = R"(

    #version 330 core 

    out vec3 color;

    void main()
    {
	color = vec3(0.0098, 0.0098, 0.439);
    }
    )";


    template<>
    BufferedShapeDrawer<Curve>::BufferedShapeDrawer()
    : vao(), shader(curve_vert, curve_frag, nullptr) {

	vertex_attributes.push_back(0);
	vao.bind();
	vao.buffers.emplace("vert_data", GL_ARRAY_BUFFER);

    }

    template<>
    void BufferedShapeDrawer<Curve>::buffer_data(){
	vao.bind();
	vao["vert_data"].bind();
	vao["vert_data"].buffer_data(host_buffer.pts.size(), host_buffer.pts.data(), GL_STATIC_DRAW);
	vao["vert_data"].set_attr_pointer<GLfloat>(0, 3, 3, 0);
    }

    template<>
    void BufferedShapeDrawer<Curve>::draw_arrays(){
	for(int icurve = 0; icurve < host_buffer.iterators.size() - 1; ++icurve){
	    GLsizei count = host_buffer.iterators[icurve + 1] - host_buffer.iterators[icurve];
	    glDrawArrays(GL_LINE_STRIP, host_buffer.iterators[icurve], count);
	}
    }
}
