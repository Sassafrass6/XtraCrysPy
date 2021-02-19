from dataclasses import dataclass

@dataclass
class Link:
    node: str
    target_index: int

material_data = {
    'mat.glossy': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.opaque': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.metallic': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.4, 0.0, 0.0, 0.3, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.glass': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.1, 1.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.translucent': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.1, 0.9, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.plastic1': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.1, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.plastic2': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, None, None, 0.0, 0.1, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, None, None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        }
    },
    'mat.chalk': {
        'OutputMaterial.0': {
            'inputs': (
                Link('BsdfPrincipled.0', 0), None, None,
            )
        },
        'BsdfPrincipled.0': {
            'inputs': (
                None, None, (1.0, 0.2, 0.1), None, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0,
                0.0, 0.5, 0.0, 0.03, 1.45, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0), 1.0, Link('Bump.0', 0), None, None,
            ),
            'outputs': (
                Link('OutputMaterial.0', 0),
            )
        },
        'Bump.0': {
            'invert': True,
            'inputs': (
                0.2, 1.0, Link('TexNoise.0', 1),
            ),
            'outputs': (
                Link('BsdfPrincipled.0', 20),
            )
        },
        'TexNoise.0': {
            'noise_dimensions': '3D',
            'inputs': (
                None, None, 1.0, 20.0, 0.9, 0.0,
            ),
            'outputs': (
                None, Link('Bump.0', 2),
            )
        }
    },
    
}