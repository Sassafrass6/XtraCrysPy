from fury.lib import Actor2D,PolyDataMapper2D
from fury import ui

import vtkmodules.vtkFiltersSources as fsvtk
LineSource = fsvtk.vtkLineSource

# Create my Line2D class, until I'm allowed to merge the branch into fury
class Line2D(ui.core.UI):
    """A 2D line UI component."""

    def __init__(self, p1, p2, color=(1, 0, 0), opacity=1.0, width=2):
        """Initialize a 2D Line.

        Parameters
        ----------
        p1 : (int)
            First endpoint of the line
        p2 : (int)
            Second endpoint of the line
        color : (float, float, float), optional
            Must take values in [0, 1].
        opacity : float, optional
            Must take values in [0, 1].
        width : float
            Thickness of the line

        """
        super(Line2D, self).__init__()

        self.color = color
        self.lwidth = width
        self.opacity = opacity

        self.lp1 = p1
        self.lp2 = p2

        self._line.SetPoint1((*p1, 0))
        self._line.SetPoint2((*p2, 0))
        self._line.Update()

        self.actor.GetProperty().SetColor(*color)
        self.actor.GetProperty().SetLineWidth(width)

    def _setup(self):
        """Setup this UI component.

        Creating the line actor used internally.

        """
        # Setting up line actor.
        self._line = LineSource()

        # Mapper
        mapper = PolyDataMapper2D()
        mapper.SetInputConnection(self._line.GetOutputPort())

        # Actor
        self.actor = Actor2D()
        self.actor.SetMapper(mapper)

        # Add default events listener to the VTK actor.
        self.handle_events(self.actor)

    def _get_actors(self):
        """Get the actors composing this UI component."""
        return [self.actor]

    def _add_to_scene(self, scene):
        """Add all subcomponents or VTK props that compose this UI component.

        Parameters
        ----------
        scene : scene

        """
        scene.add(self.actor)

    def _get_size(self):
        size = np.abs(self.lp2 - self.lp1)
        return size

    def _set_position(self, coords):
        """Set the lower-left corner position of this UI bounding box.
           Line requires two points, so one coordinate is useless

        Parameters
        ----------
        coords: (float, float)
            Absolute pixel coordinates (x, y).
        """
        pass

    @property
    def color(self):
        """Get the color of this UI component."""
        color = self.actor.GetProperty().GetColor()
        return np.asarray(color)

    @color.setter
    def color(self, color):
        """Set the color of this UI component.

        Parameters
        ----------
        color : (float, float, float)
            RGB. Must take values in [0, 1].

        """
        self.actor.GetProperty().SetColor(*color)

    @property
    def opacity(self):
        """Get the opacity of this UI component."""
        return self.actor.GetProperty().GetOpacity()

    @opacity.setter
    def opacity(self, opacity):
        """Set the opacity of this UI component.

        Parameters
        ----------
        opacity : float
            Degree of transparency. Must be between [0, 1].

        """
        self.actor.GetProperty().SetOpacity(opacity)

    @property
    def p1(self):
      return self._line.GetPoint1()[:-1]

    @p1.setter
    def p1(self, p1):
      self._line.SetPoint1((*p1, 0))
      self._line.Update()

    @property
    def p2(self):
      return self._line.GetPoint2()[:-1]

    @p2.setter
    def p2(self, p2):
      self._line.SetPoint2((*p2, 0))
      self._line.Update()

    @property
    def width(self):
      return self._line.GetLineWidth()

    @width.setter
    def width(self, width):
      self.actor.GetProperty().SetLineWidth(width)

